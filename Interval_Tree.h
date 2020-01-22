#ifndef INTERVAL_TREE_H
#define INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>
#include <cassert>

template <class Scalar, typename Value>
class Interval 
{
    public:
        Scalar start;
        Scalar stop;
        Value value;

        Interval(const Scalar& s, const Scalar& e, const Value& v) : start(std::min(s, e)), stop(std::max(s, e)), value(v) {}
};

template <class Scalar, typename Value>
Value intervalStart(const Interval<Scalar, Value>& i) { return i.start; }

template <class Scalar, typename Value>
Value intervalStop(const Interval<Scalar, Value>& i) { return i.stop; }

template <class Scalar, typename Value>
std::ostream& operator<<(std::ostream& out, const Interval<Scalar, Value>& i)
{
    out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class Scalar, typename Value>
class Interval_Tree 
{
    public:
        typedef Interval<Scalar, Value> interval;
        typedef std::vector<interval> interval_vector;

        struct Interval_Start_Cmp
        {
            bool operator()(const interval& a, const interval& b) { return a.start < b.start; }
        };

        struct Interval_Stop_Cmp
        {
            bool operator()(const interval& a, const interval& b) { return a.stop < b.stop; }
        };

        Interval_Tree() : left(nullptr), right(nullptr), center(0) {}
        ~Interval_Tree() = default;

        std::unique_ptr<Interval_Tree> clone() const { return std::unique_ptr<Interval_Tree>(new Interval_Tree(*this)); }

        Interval_Tree(const Interval_Tree& other)
            : intervals(other.intervals),
              left(other.left ? other.left->clone() : nullptr),
              right(other.right ? other.right->clone() : nullptr),
              center(other.center) {}

        Interval_Tree& operator=(Interval_Tree&&) = default;
        Interval_Tree(Interval_Tree&&) = default;

        Interval_Tree& operator=(const Interval_Tree& other)
        {
            center = other.center;
            intervals = other.intervals;
            left = other.left ? other.left->clone() : nullptr;
            right = other.right ? other.right->clone() : nullptr;

            return *this;
        }

        Interval_Tree(
            interval_vector&& ivals,
            std::size_t depth = 16,
            std::size_t minBucket = 64,
            std::size_t maxBucket = 512,
            Scalar leftExtent = 0,
            Scalar rightExtent = 0)
                : left(nullptr),
                  right(nullptr)
        {
            depth--;

            const auto minimaxStart = std::minimax_element(ivals.begin(), ivals.end(), Interval_Start_Cmp());
            const auto minimaxStop = std::minimax_element(ivals.begin(), ivals.end(), Interval_Stop_Cmp());

            if (!ivals.empty())
            {
                center = (minimaxStart.first->start + minimaxStop.second->stop) / 2;
            }

            if (leftExtent == 0 && rightExtent == 0)
            {
                // Sort intervals by start
                std::sort(ivals.begin(), ivals.end(), Interval_Start_Cmp());
            }
            else 
            {
                assert(std::is_sorted(ivals.begin(), ivals.end(), Interval_Start_Cmp()));
            }

            if (depth == 0 || (ivals.size() < minBucket && ivals.size() < maxBucket))
            {
                std::sort(ivals.begin(), ivals.end(), Interval_Start_Cmp());
                intervals = std::move(ivals);
                assert(is_valid().first);

                return;
            }
            else 
            {
                Scalar leftp = 0;
                Scalar rightp = 0;

                if (leftExtent || rightExtent)
                {
                    leftp = leftExtent;
                    rightp = rightExtent;
                }
                else 
                {
                    leftp = ivals.front().start;
                    rightp = std::max_element(ivals.begin(), ivals.end(), Interval_Stop_Cmp())->stop;
                }

                interval_vector lefts;
                interval_vector rights;

                for (typename interval_vector::const_iterator i = ivals.begin(); i != ivals.end(); i++)
                {
                    const interval& interval = *i;

                    if (interval.stop < center)
                    {
                        lefts.push_back(interval);
                    }
                    else if (interval.start > center)
                    {
                        rights.push_back(interval);
                    }
                    else 
                    {
                        assert(interval.start <= center);
                        assert(center <= interval.stop);
                        intervals.push_back(interval);
                    }
                }

                if (!lefts.empty())
                {
                    left.reset(new Interval_Tree(std::move(lefts), depth, minBucket, maxBucket, leftp, center));
                }

                if (!rights.empty())
                {
                    right.reset(new Interval_Tree(std::move(rights), depth, minBucket, maxBucket, center, rightp));
                }
            }

            assert(is_valid().first);
        }

        // Call f on all intervals near the range [start, stop]
        template <class UnaryFunction>
        void visit_near(const Scalar& start, const Scalar& stop, UnaryFunction f) const
        {
            if (!intervals.empty() && ! (stop < intervals.front().start))
            {
                for (auto & i : intervals)
                {
                    f(i);
                }
            }

            if (left && start <= center)
            {
                left->visit_near(start, stop, f);
            }

            if (right && stop >= center)
            {
                right->visit_near(start, stop, f);
            }
        }

        // Call f on all intervals crossing pos
        template <class UnaryFunction>
        void visit_crossing(const Scalar& pos, UnaryFunction f) const
        {
            visit_crossing(pos, pos, f);
        }

        // Call f on all intervals overlapping [start, stop]
        template <class UnaryFunction>
        void visit_overlapping(const Scalar& start, const Scalar& stop, UnaryFunction f) const
        {
            auto filterF = [&](const interval& interval)
            {
                if (interval.stop >= start && interval.start <= stop)
                {
                    // Only apply f if overlapping
                    f(interval);
                }
            };

            visit_near(start, stop, filterF);
        }

        // Call f on all intervals contained within [start, stop]
        template <class UnaryFunction>
        void visit_contained(const Scalar& start, const Scalar& stop, UnaryFunction f) const
        {
            auto filterF = [&](const interval& interval)
            {
                if (start <= interval.start && interval.stop <= stop)
                {
                    f(interval);
                }
            };

            visit_near(start, stop, filterF);
        }

        interval_vector findOverlapping(const Scalar& start, const Scalar& stop) const
        {
            interval_vector result;

            visit_overlapping(start, stop,
                              [&](const interval& interval)
                              {
                                  result.emplace_back(interval);
                              });
            
            return result;
        }

        interval_vector findContained(const Scalar& start, const Scalar& stop) const
        {
            interval_vector result;

            visit_contained(start, stop,
                            [&](const interval& interval)
                            {
                                result.push_back(interval);
                            });

            return result;
        }

        bool empty() const
        {
            if (left && !left->empty())
                return false;

            if (!intervals.empty())
                return false;

            if (right && !right->empty())
                return false;

            return true;
        }

        template <class UnaryFunction>
        void visit_all(UnaryFunction f) const
        {
            if (left)
            {
                left->visit_all(f);
            }

            std::for_each(intervals.begin(), intervals.end(), f);

            if (right)
            {
                right->visit_all(f);
            }
        }

        std::pair<Scalar, Scalar> extentBruteForce() const
        {
            struct Extent 
            {
                std::pair<Scalar, Scalar> x = 
                { 
                    std::numeric_limits<Scalar>::max(), 
                    std::numeric_limits<Scalar>::min() 
                };

                void operator()(const interval & interval)
                {
                    x.first = std::min(x.first, interval.start);
                    x.second = std::max(x.second, interval.stop);
                }
            };  

            Extent extent;

            visit_all([&](const interval & interval) 
            {
                extent(interval);
            });

            return extent.x;
        }

        // Check all constraints
        // If first is false, second is invalid
        std::pair<bool, std::pair<Scalar, Scalar>> is_valid() const
        {
            const auto minimaxStart = std::minimax_element(intervals.begin(), intervals.end(), Interval_Start_Cmp());
            const auto minimaxStop = std::minimax_element(intervals.begin(), intervals.end(), Interval_Stop_Cmp());

            std::pair<bool, std::pair<Scalar, Scalar> result =
            {
                std::numeric_limits<Scalar>::max();
                std::numeric_limits<Scalar>::min();
            };

            if (!intervals.empty())
            {
                result.second.first = std::min(result.second.first, minimaxStart.first->start);
                result.second.second = std::min(result.second.second, minimaxStop.second->stop);
            }

            if (left)
            {
                auto valid = left->is_valid();
                result.first & = valid.first;

                result.second.first = std::min(result.second.first, valid.second.first);
                result.second.second = std::min(result.second.second, valid.second.second);

                if (!result.filterF)
                    return result;

                if (valid.second.first <= center)
                {
                    result.first = false;
                    return result;
                }
            }

            if (right)
            {
                auto valid = right->is_valid();
                result.first & = valid.first;

                result.second.first = std::min(result.second.first, valid.second.first);
                result.second.second = std::min(result.second.second, valid.second.second);

                if (!result.first)
                    return result;

                if (valid.second.first <= center)
                {
                    result.first = false;
                    return result;
                }
            }

            if (!std::is_sorted(intervals.begin(), intervals.end(), Interval_Start_Cmp()))
            {
                result.first = false;
            }

            return result;
        }

        friend std::ostream& writeOut(std::ostream& os, const Interval_Tree& itree, std::size_t depth = 0)
        {
            auto pad = [&]()
            {
                for (std::size_t i = 0; i != depth; i++)
                {
                    os << ' ';
                }
            };

            pad();

            os << "center: " << itree.center << '\n';

            if (itree.left)
            {
                pad();

                os << "left:\n";

                writeOut(os, *itree.left, depth + 1);
            }
            else 
            {
                pad();

                os << "left: nullptr\n";
            }

            if (itree.right)
            {
                pad();

                os << "right:\n";

                writeOut(os, *itree.right, depth + 1);
            }
            else 
            {
                pad();

                os << "right: nullptr\n";
            }

            return os;
        }

    private:
        interval_vector intervals;
        std::unique_ptr<Interval_Tree> left;
        std::unique_ptr<Interval_Tree> right;
        Scalar center;
};

#endif /* INTERVAL_TREE_H */