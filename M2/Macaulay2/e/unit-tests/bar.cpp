//
// Created by Radoslav on 2019-02-21.
//

#include <gtest/gtest.h>

static constexpr bool bar_assert(bool flag) {
    return flag;
}

TEST(Test02, uniquetest) {
    EXPECT_TRUE(bar_assert(true));
}
