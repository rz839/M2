//
// Created by Radoslav on 2019-02-21.
//

#include "foo.h"
#include <gtest/gtest.h>

TEST(Test01, subtest01)
{
    EXPECT_EQ(scratch::foo(5), 25);
}

TEST(Test01, subtest02) {
    EXPECT_GE(scratch::foo(10),99);
}

