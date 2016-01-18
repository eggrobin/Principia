﻿#pragma once

namespace principia {
namespace base {

// This class helps in declaring that a type is "mappable", i.e., that the maps
// declared in principia::geometry can be act on it through the operator().  To
// use it, declare a specialization for the proper |Functor| (a class that must
// have an operator()) and |T| (the class to be made mappable).  The third
// parameter is an enabler and maybe used to restrict specializations with
// SFINAE.
// Any specialization must export two declarations in its public part:
//   o A type named |type|, which is the result of applying |Functor| to |T|.
//   o A static function named |Do|, which applies a functor to a value of |T|
//     and returns a |type|.
// See the comments below for example declarations.
template<typename Functor, typename T, typename = void>
class Mappable {
 public:
  // using type = void;
  // static type Do(Functor const& functor, T const& t);
};

}  // namespace base
}  // namespace principia
