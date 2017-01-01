cc_library(
    name = "macros",
    hdrs = ["macros.hpp"],
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
