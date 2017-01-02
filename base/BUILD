cc_library(
    name = "array",
    hdrs = ["array.hpp"],
    srcs = ["array_body.hpp"],
)

cc_library(
    name = "macros",
    hdrs = ["macros.hpp"],
)

cc_library(
    name = "not_constructible",
    hdrs = ["not_constructible.hpp"],
)

cc_library(
    name = "not_null",
    hdrs = ["not_null.hpp"],
    srcs = ["not_null_body.hpp"],
    deps = [
        ":not_constructible",
    ],
)

cc_test(
    name = "not_null_test",
    srcs = ["not_null_test.cpp"],
    deps = [
        ":not_null",
        "//deps:gtest",
        "//deps:gmock_main",
    ],
)

cc_library(
    name = "pull_serializer",
    hdrs = ["pull_serializer.hpp"],
    srcs = ["pull_serializer_body.hpp"],
    deps = [
        ":array",
        ":not_null",
    ],
)

cc_test(
    name = "pull_serializer_test",
    srcs = ["pull_serializer_test.cpp"],
    deps = [
        ":macros",
        ":pull_serializer",
        "//deps:gtest",
        "//deps:gmock_main",
        "//deps/protobuf:protobuf",
    ],
)

cc_library(
    name = "status",
    srcs = ["status.cpp"],
    hdrs = ["status.hpp"],
    deps = [
        ":macros",
        "//deps:glog",
    ],
)

cc_test(
    name = "status_test",
    srcs = ["status_test.cpp"],
    deps = [
        ":status",
        "//deps:gtest",
        "//deps:gmock_main",
    ],
)
