/**
 * @file qm_worker_test.cc
 */

#include <gtest/gtest.h>

#include "../stdheader.h"

#include "qm_worker_test.h"

#include "../interaction/qmmm/qmmm_interaction.h"

namespace testing {

QM_Worker_Test::QM_Worker_Test(const Parameter& parameter) : test_sim_(parameter) {}

}