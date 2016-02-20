﻿
#pragma once

#include "physics/discrete_trajectory.hpp"

#include <algorithm>
#include <list>
#include <map>

#include "geometry/named_quantities.hpp"
#include "glog/logging.h"

namespace principia {

using base::make_not_null_unique;
using geometry::Instant;

namespace physics {
namespace internal {

template<typename Frame>
Instant const& ForkableTraits<DiscreteTrajectory<Frame>>::time(
    TimelineConstIterator const it) {
  return it->first;
}

template<typename Frame>
Instant const& DiscreteTrajectoryIterator<Frame>::time() const {
  return this->current()->first;
}

template<typename Frame>
DegreesOfFreedom<Frame> const&
DiscreteTrajectoryIterator<Frame>::degrees_of_freedom() const {
  return this->current()->second;
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame>*>
DiscreteTrajectoryIterator<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectoryIterator<Frame> const*>
DiscreteTrajectoryIterator<Frame>::that() const {
  return this;
}

}  // namespace internal

template<typename Frame>
DiscreteTrajectory<Frame>::~DiscreteTrajectory() {
  if (on_destroy_) {
    on_destroy_(this);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::set_on_destroy(
    OnDestroyCallback on_destroy) {
  on_destroy_ = on_destroy;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::OnDestroyCallback
DiscreteTrajectory<Frame>::get_on_destroy() const {
  return on_destroy_;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::Iterator
DiscreteTrajectory<Frame>::last() const {
  return --this->End();
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkWithCopy(Instant const& time) {
  // May be at |timeline_end()| if |time| is the fork time of this object.
  auto timeline_it = timeline_.find(time);
  CHECK(timeline_it != timeline_end() ||
        (!this->is_root() && time == this->Fork().time()))
      << "NewFork at nonexistent time " << time;

  auto const fork = this->NewFork(timeline_it);

  // Copy the tail of the trajectory in the child object.
  if (timeline_it != timeline_.end()) {
    fork->timeline_.insert(++timeline_it, timeline_.end());
  }
  return fork;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*>
DiscreteTrajectory<Frame>::NewForkAtLast() {
  auto end = timeline_.end();
  if (timeline_.empty()) {
    return this->NewFork(end);
  } else {
    return this->NewFork(--end);
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::Append(
    Instant const& time,
    DegreesOfFreedom<Frame> const& degrees_of_freedom) {
  CHECK(this->is_root() || time > this->Fork().time())
       << "Append at " << time << " which is before fork time "
       << this->Fork().time();

  if (!timeline_.empty() && timeline_.cbegin()->first == time) {
    LOG(WARNING) << "Append at existing time " << time
                 << ", time range = [" << this->Begin().time() << ", "
                 << last().time() << "]";
    return;
  }
  auto it = timeline_.emplace_hint(timeline_.end(),
                                   time,
                                   degrees_of_freedom);
  // Decrementing |end()| is much faster than incrementing |it|.  Don't ask.
  CHECK(--timeline_.end() == it) << "Append out of order";
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetAfter(Instant const& time) {
  this->DeleteAllForksAfter(time);

  // Get an iterator denoting the first entry with time > |time|.  Remove that
  // entry and all the entries that follow it.  This preserves any entry with
  // time == |time|.
  auto const it = timeline_.upper_bound(time);
  timeline_.erase(it, timeline_.end());
}

template<typename Frame>
void DiscreteTrajectory<Frame>::ForgetBefore(Instant const& time) {
  this->DeleteAllForksBefore(time);

  // Get an iterator denoting the first entry with time >= |time|.  Remove all
  // the entries that precede it.  This preserves any entry with time == |time|.
  auto it = timeline_.lower_bound(time);
  timeline_.erase(timeline_.begin(), it);
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteToMessage(
    not_null<serialization::Trajectory*> const message) const {
  LOG(INFO) << __FUNCTION__;
  CHECK(this->is_root());
  WriteSubTreeToMessage(message);
  LOG(INFO) << NAMED(this);
  LOG(INFO) << NAMED(message->SpaceUsed());
  LOG(INFO) << NAMED(message->ByteSize());
}

template<typename Frame>
not_null<std::unique_ptr<DiscreteTrajectory<Frame>>>
DiscreteTrajectory<Frame>::ReadFromMessage(
    serialization::Trajectory const& message) {
  auto trajectory = make_not_null_unique<DiscreteTrajectory>();
  trajectory->FillSubTreeFromMessage(message);
  return trajectory;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame>*> DiscreteTrajectory<Frame>::that() {
  return this;
}

template<typename Frame>
not_null<DiscreteTrajectory<Frame> const*>
DiscreteTrajectory<Frame>::that() const {
  return this;
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_begin() const {
  return timeline_.begin();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_end() const {
  return timeline_.end();
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_find(Instant const& time) const {
  return timeline_.find(time);
}

template<typename Frame>
typename DiscreteTrajectory<Frame>::TimelineConstIterator
DiscreteTrajectory<Frame>::timeline_lower_bound(Instant const& time) const {
  return timeline_.lower_bound(time);
}

template<typename Frame>
bool DiscreteTrajectory<Frame>::timeline_empty() const {
  return timeline_.empty();
}

template<typename Frame>
void DiscreteTrajectory<Frame>::WriteSubTreeToMessage(
    not_null<serialization::Trajectory*> const message) const {
  Forkable<DiscreteTrajectory, Iterator>::WriteSubTreeToMessage(message);
  for (auto const& pair : timeline_) {
    Instant const& instant = pair.first;
    DegreesOfFreedom<Frame> const& degrees_of_freedom = pair.second;
    auto const instantaneous_degrees_of_freedom = message->add_timeline();
    instant.WriteToMessage(instantaneous_degrees_of_freedom->mutable_instant());
    degrees_of_freedom.WriteToMessage(
        instantaneous_degrees_of_freedom->mutable_degrees_of_freedom());
  }
}

template<typename Frame>
void DiscreteTrajectory<Frame>::FillSubTreeFromMessage(
    serialization::Trajectory const& message) {
  for (auto timeline_it = message.timeline().begin();
       timeline_it != message.timeline().end();
       ++timeline_it) {
    Append(Instant::ReadFromMessage(timeline_it->instant()),
           DegreesOfFreedom<Frame>::ReadFromMessage(
               timeline_it->degrees_of_freedom()));
  }
  Forkable<DiscreteTrajectory, Iterator>::FillSubTreeFromMessage(message);
}

}  // namespace physics
}  // namespace principia
