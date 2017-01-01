cc_library(
    name = "glog",
    srcs = ["glog/.libs/libglog.a"],
    hdrs = glob(["glog/src/glog/*.h"]),
    includes = ["glog/src"],
    visibility = ["//visibility:public"],
    linkopts = ["-pthread"],
)

cc_library(
    name = "gtest",
    srcs = glob(
        [
            "googletest/src/*.cc",
            "googletest/src/*.h",
        ],
        exclude = [
            "googletest/src/gtest-all.cc",
            "googletest/src/gtest_main.cc",
        ]
    ),
    hdrs = glob(["googletest/include/**/*.h"]),
    includes = ["googletest/include"],
    visibility = ["//visibility:public"],
    copts = ["-Ideps/googletest"]
)

cc_library(
    name = "gmock",
    srcs = [
        "googletest/src/gmock-all.cc",
    ],
    hdrs = glob(["googlemock/include/**/*.h"]),
    deps = [":gtest"],
    includes = ["googlemock/include"],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "gmock_main",
    srcs = ["googlemock/src/gmock_main.cc"],
    deps = [
        ":glog",
        ":gmock"
    ],
)
