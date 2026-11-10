/*  This file is part of the OpenLB library
*
*  Copyright (C) 2026 Adrian Kummerlaender
*  E-mail contact: info@openlb.net
*  The most recent release of OpenLB can be downloaded at
*  <http://www.openlb.net/>
*
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License
*  as published by the Free Software Foundation; either version 2
*  of the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public
*  License along with this program; if not, write to the Free
*  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
*  Boston, MA  02110-1301, USA.
*/

#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace olb;

namespace olb::parameters {

struct TEST_SCALAR   : public descriptors::FIELD_BASE<1> { };
struct TEST_VECTOR   : public descriptors::FIELD_BASE<0,1> { };
struct TEST_BASE     : public descriptors::FIELD_BASE<1> { };
struct TEST_DERIVED  : public descriptors::FIELD_BASE<1> { };
struct TEST_OTHER    : public descriptors::FIELD_BASE<1> { };

struct TEST_COUNT : public descriptors::TYPED_FIELD_BASE<int,1> { };

struct TEST_POSITIVE : public descriptors::FIELD_BASE<1> {
  template <typename T, typename DESCRIPTOR, typename FIELD>
  static constexpr auto isValid(FieldD<T,DESCRIPTOR,FIELD> value) {
    return value > 0;
  }
};

/// Only ever read through the overlay of unrelated fields, to force map growth
struct TEST_FILLER_0 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_1 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_2 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_3 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_4 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_5 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_6 : public descriptors::FIELD_BASE<1> { };
struct TEST_FILLER_7 : public descriptors::FIELD_BASE<1> { };

}

using T = double;
using DESCRIPTOR = descriptors::D2Q9<>;
using Parameters = ParametersD<T,DESCRIPTOR>;

namespace {

/// Replace the process-global override sources by exactly the given CLI tokens.
/// A --case-file pointing nowhere keeps a stray case.xml out of the test.
void setOverrides(std::vector<std::string> tokens)
{
  tokens.insert(tokens.begin(), "testParametersD");
  std::vector<char*> argv;
  for (auto& token : tokens) {
    argv.push_back(token.data());
  }
  int argc = argv.size();
  singleton::input().load(argc, argv.data());
}

void clearOverrides()
{
  setOverrides({"--case-file", "/nonexistent/case.xml"});
}

struct OverrideGuard {
  OverrideGuard() { clearOverrides(); }
  ~OverrideGuard() { clearOverrides(); }
};

}


TEST(parametersD, eagerScalarAndVector)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_SCALAR>(0.25);
  parameters.set<parameters::TEST_VECTOR>(Vector<T,2>{1.0, 2.0});

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 0.25, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_VECTOR>()[0], 1.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_VECTOR>()[1], 2.0, 1e-12);
}

TEST(parametersD, providesAndFallback)
{
  OverrideGuard guard;
  Parameters parameters;

  EXPECT_FALSE(parameters.provides<parameters::TEST_SCALAR>());
  EXPECT_NEAR(parameters.getOrFallback<parameters::TEST_SCALAR>(7.0), 7.0, 1e-12);

  parameters.set<parameters::TEST_SCALAR>(0.5);

  EXPECT_TRUE(parameters.provides<parameters::TEST_SCALAR>());
  EXPECT_NEAR(parameters.getOrFallback<parameters::TEST_SCALAR>(7.0), 0.5, 1e-12);
}

TEST(parametersD, versionBumpsOnSet)
{
  OverrideGuard guard;
  Parameters parameters;

  const std::size_t initial = parameters.version();
  parameters.set<parameters::TEST_SCALAR>(1.0);
  EXPECT_EQ(parameters.version(), initial + 1);
  parameters.set<parameters::TEST_SCALAR>(2.0);
  EXPECT_EQ(parameters.version(), initial + 2);
}

TEST(parametersD, dependentTracksItsBase)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(1.0);
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 2.0 * p.get<parameters::TEST_BASE>();
  });

  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 2.0, 1e-12);

  parameters.set<parameters::TEST_BASE>(3.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 6.0, 1e-12);
}

TEST(parametersD, setOverridesDependency)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(1.0);
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 2.0 * p.get<parameters::TEST_BASE>();
  });
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 2.0, 1e-12);

  // A manual set drops the dependency for good
  parameters.set<parameters::TEST_DERIVED>(99.0);
  parameters.set<parameters::TEST_BASE>(5.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 99.0, 1e-12);
}

/// Regression: a caught cyclic-dependency error used to leave the field flagged
/// as under computation, so every later access threw the same error forever
TEST(parametersD, cyclicDependencyThrowsAndRecovers)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>([](Parameters& p) {
    return p.get<parameters::TEST_DERIVED>() + 1.0;
  });
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return p.get<parameters::TEST_BASE>() + 1.0;
  });

  EXPECT_THROW(parameters.get<parameters::TEST_BASE>(), std::exception);

  // Breaking the cycle must make both parameters usable again
  parameters.set<parameters::TEST_DERIVED>(1.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_BASE>(), 2.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 1.0, 1e-12);
}

/// The entry reference held across calculate() must survive the insertions that
/// evaluating the lambda performs
TEST(parametersD, entryReferenceSurvivesMapGrowth)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    T sum = 0;
    sum += p.get<parameters::TEST_FILLER_0>();
    sum += p.get<parameters::TEST_FILLER_1>();
    sum += p.get<parameters::TEST_FILLER_2>();
    sum += p.get<parameters::TEST_FILLER_3>();
    sum += p.get<parameters::TEST_FILLER_4>();
    sum += p.get<parameters::TEST_FILLER_5>();
    sum += p.get<parameters::TEST_FILLER_6>();
    sum += p.get<parameters::TEST_FILLER_7>();
    return sum + 42.0;
  });

  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 42.0, 1e-12);
}

TEST(parametersD, subgroupsAreIndependent)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters["reporter"].set<parameters::TEST_SCALAR>(2.0);
  parameters["reporter"]["vtk"].set<parameters::TEST_SCALAR>(3.0);

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
  EXPECT_NEAR(parameters["reporter"].get<parameters::TEST_SCALAR>(), 2.0, 1e-12);
  EXPECT_NEAR(parameters["reporter"]["vtk"].get<parameters::TEST_SCALAR>(), 3.0, 1e-12);
}

TEST(parametersD, subgroupReachesRoot)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(4.0);
  parameters["reporter"]["timer"].set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 2.0 * p.root().get<parameters::TEST_BASE>();
  });

  EXPECT_NEAR(parameters["reporter"]["timer"].get<parameters::TEST_DERIVED>(), 8.0, 1e-12);
  EXPECT_NEAR(parameters["reporter"].parent().get<parameters::TEST_BASE>(), 4.0, 1e-12);
}

/// A copy must be rooted in itself, never in the tree it came from
TEST(parametersD, copiedTreeIsRootedInItself)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(4.0);
  parameters["reporter"]["timer"].set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 2.0 * p.root().get<parameters::TEST_BASE>();
  });

  Parameters copy = parameters;
  copy.set<parameters::TEST_BASE>(10.0);

  EXPECT_NEAR(copy["reporter"]["timer"].get<parameters::TEST_DERIVED>(), 20.0, 1e-12);
  EXPECT_NEAR(parameters["reporter"]["timer"].get<parameters::TEST_DERIVED>(), 8.0, 1e-12);
}

TEST(parametersD, movedTreeIsRootedInItself)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(4.0);
  parameters["reporter"]["timer"].set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 2.0 * p.root().get<parameters::TEST_BASE>();
  });

  Parameters moved = std::move(parameters);
  EXPECT_NEAR(moved["reporter"]["timer"].get<parameters::TEST_DERIVED>(), 8.0, 1e-12);
  EXPECT_EQ(&moved["reporter"].parent(), &moved);
}

/// A subgroup copied out of its tree has no root, and must say so
TEST(parametersD, detachedSubgroupHasNoRoot)
{
  OverrideGuard guard;
  Parameters parameters;
  parameters.set<parameters::TEST_BASE>(4.0);

  Parameters detached = parameters["reporter"];

  EXPECT_THROW(detached.root(), std::runtime_error);
  EXPECT_THROW(detached.parent(), std::runtime_error);
  EXPECT_EQ(&parameters.root(), &parameters);
}

TEST(parametersD, validationRejectsInvalidValue)
{
  OverrideGuard guard;
  Parameters parameters;

  EXPECT_THROW(parameters.set<parameters::TEST_POSITIVE>(-1.0), std::invalid_argument);

  parameters.set<parameters::TEST_POSITIVE>([](Parameters&) { return -2.0; });
  EXPECT_THROW(parameters.get<parameters::TEST_POSITIVE>(), std::invalid_argument);

  // The failed dependent evaluation must not poison the parameter
  parameters.set<parameters::TEST_POSITIVE>(4.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_POSITIVE>(), 4.0, 1e-12);
}


/// The rows are only safe to copy because they cannot bind to an instance.
/// A capturing lambda must not be convertible to the row type.
TEST(parametersD, capturingLambdaIsNotAValidRow)
{
  using Row = meta::function<void(int)>;

  int captured = 0;
  auto capturing = [captured](int) { (void) captured; };
  auto captureless = [](int) { };

  static_assert(!std::is_constructible_v<Row, decltype(capturing)>,
                "a capturing lambda must not be storable in a type-erased row");
  static_assert(std::is_constructible_v<Row, decltype(captureless)>,
                "a captureless lambda must remain storable");

  EXPECT_FALSE((std::is_constructible_v<Row, decltype(capturing)>));
  EXPECT_TRUE((std::is_constructible_v<Row, decltype(captureless)>));
}

/// Derivations are rejected at compile time unless they capture nothing,
/// which is what lets a copy or a conversion evaluate them correctly
TEST(parametersD, capturingDerivationIsNotStorable)
{
  int captured = 0;
  auto capturing = [captured](Parameters&) { return 1.0 * captured; };
  auto captureless = [](Parameters&) { return 1.0; };

  static_assert(!std::is_empty_v<decltype(capturing)>,
                "a capturing lambda must not satisfy the derivation requirement");
  static_assert(std::is_empty_v<decltype(captureless)>,
                "a captureless lambda must satisfy the derivation requirement");

  EXPECT_FALSE(std::is_empty_v<decltype(capturing)>);
  EXPECT_TRUE(std::is_empty_v<decltype(captureless)>);
}

TEST(parametersD, isCopyableAndMovable)
{
  static_assert(std::is_copy_constructible_v<Parameters>);
  static_assert(std::is_copy_assignable_v<Parameters>);
  static_assert(std::is_move_constructible_v<Parameters>);
  static_assert(std::is_move_assignable_v<Parameters>);

  EXPECT_TRUE(std::is_copy_constructible_v<Parameters>);
  EXPECT_TRUE(std::is_move_constructible_v<Parameters>);
}

TEST(parametersD, unparseableCliOverrideThrows)
{
  OverrideGuard guard;
  setOverrides({"--TEST_SCALAR", "not-a-number"});

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  EXPECT_THROW(parameters.get<parameters::TEST_SCALAR>(), std::invalid_argument);
}

TEST(parametersD, unparseableXmlOverrideThrows)
{
  OverrideGuard guard;

  const std::string filename = "testParametersD_bad.xml";
  {
    std::ofstream out(filename, std::ios::trunc);
    out << "<?xml version=\"1.0\"?>\n"
        << "<Case>\n<Parameters>\n"
        << "  <TEST_SCALAR>not-a-number</TEST_SCALAR>\n"
        << "</Parameters>\n</Case>\n";
  }
  setOverrides({"--case-file", filename});

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  EXPECT_THROW(parameters.get<parameters::TEST_SCALAR>(), std::invalid_argument);

  std::remove(filename.c_str());
}

/// An override may not smuggle a value past the field's own validation
TEST(parametersD, overrideStillValidates)
{
  OverrideGuard guard;
  setOverrides({"--TEST_POSITIVE", "-3.0"});

  Parameters parameters;
  parameters.set<parameters::TEST_POSITIVE>(1.0);

  EXPECT_THROW(parameters.get<parameters::TEST_POSITIVE>(), std::invalid_argument);
}

/// The recursion guard must be cleared on the way out, so a second access
/// reports the lambda's own error again instead of a phantom cycle
TEST(parametersD, failingDependencyReportsSameErrorTwice)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(
    [](Parameters&) -> FieldD<T,DESCRIPTOR,parameters::TEST_BASE> {
      throw std::runtime_error("deliberate");
    });

  auto messageOf = [&parameters]() {
    try {
      parameters.get<parameters::TEST_BASE>();
    } catch (const std::exception& e) {
      return std::string(e.what());
    }
    return std::string("no exception");
  };

  const std::string first = messageOf();
  const std::string second = messageOf();

  EXPECT_EQ(first, "deliberate");
  EXPECT_EQ(second, first);
}

TEST(parametersD, selfReferentialDependencyThrows)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>([](Parameters& p) {
    return p.get<parameters::TEST_BASE>() + 1.0;
  });

  EXPECT_THROW(parameters.get<parameters::TEST_BASE>(), std::runtime_error);
}

TEST(parametersD, threeWayCycleThrowsAndRecovers)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>([](Parameters& p) {
    return p.get<parameters::TEST_DERIVED>() + 1.0;
  });
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return p.get<parameters::TEST_OTHER>() + 1.0;
  });
  parameters.set<parameters::TEST_OTHER>([](Parameters& p) {
    return p.get<parameters::TEST_BASE>() + 1.0;
  });

  EXPECT_THROW(parameters.get<parameters::TEST_BASE>(), std::runtime_error);

  parameters.set<parameters::TEST_OTHER>(1.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_BASE>(), 3.0, 1e-12);
}

/// No conversion is enumerated, so the derivative count is not bounded
TEST(parametersD, copyAsArbitraryAdfDimension)
{
  OverrideGuard guard;
  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.5);

  auto many = parameters.copyAs<util::ADf<T,11>>();
  EXPECT_NEAR(many.get<parameters::TEST_SCALAR>().v(), 1.5, 1e-12);

  auto more = parameters.copyAs<util::ADf<T,32>>();
  EXPECT_NEAR(more.get<parameters::TEST_SCALAR>().v(), 1.5, 1e-12);
}


TEST(parametersD, printIsOrderedByName)
{
  OverrideGuard guard;
  Parameters parameters;

  // Deliberately set in reverse alphabetical order
  parameters.set<parameters::TEST_VECTOR>(Vector<T,2>{1.0, 2.0});
  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters.set<parameters::TEST_OTHER>(2.0);
  parameters.set<parameters::TEST_DERIVED>(3.0);
  parameters.set<parameters::TEST_BASE>(4.0);

  std::ostringstream out;
  parameters.print(out);
  const std::string text = out.str();

  const std::vector<std::string> expected{
    "TEST_BASE", "TEST_DERIVED", "TEST_OTHER", "TEST_SCALAR", "TEST_VECTOR"};

  std::size_t previous = 0;
  for (const auto& name : expected) {
    const std::size_t position = text.find(name);
    ASSERT_NE(position, std::string::npos) << name << " missing from print output";
    EXPECT_GT(position, previous) << name << " is out of order";
    previous = position;
  }
}

/// Printing evaluates dependent parameters, which may insert further entries
TEST(parametersD, printEvaluatesInsertingDependency)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    T sum = 0;
    sum += p.get<parameters::TEST_FILLER_0>();
    sum += p.get<parameters::TEST_FILLER_1>();
    sum += p.get<parameters::TEST_FILLER_2>();
    sum += p.get<parameters::TEST_FILLER_3>();
    sum += p.get<parameters::TEST_FILLER_4>();
    sum += p.get<parameters::TEST_FILLER_5>();
    sum += p.get<parameters::TEST_FILLER_6>();
    sum += p.get<parameters::TEST_FILLER_7>();
    return sum + 11.0;
  });

  std::ostringstream out;
  parameters.print(out);

  EXPECT_NE(out.str().find("TEST_DERIVED"), std::string::npos);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 11.0, 1e-12);
}

TEST(parametersD, printIsStableAcrossCopies)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_VECTOR>(Vector<T,2>{1.0, 2.0});
  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters["reporter"].set<parameters::TEST_OTHER>(2.0);

  Parameters copy = parameters;

  std::ostringstream original;
  std::ostringstream duplicate;
  parameters.print(original);
  copy.print(duplicate);

  EXPECT_EQ(original.str(), duplicate.str());
}


TEST(parametersD, copyIsIndependent)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters.set<parameters::TEST_VECTOR>(Vector<T,2>{1.0, 2.0});
  parameters["reporter"].set<parameters::TEST_SCALAR>(10.0);

  Parameters copy = parameters;

  copy.set<parameters::TEST_SCALAR>(5.0);
  copy.set<parameters::TEST_VECTOR>(Vector<T,2>{5.0, 6.0});
  copy["reporter"].set<parameters::TEST_SCALAR>(50.0);

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_VECTOR>()[0], 1.0, 1e-12);
  EXPECT_NEAR(parameters["reporter"].get<parameters::TEST_SCALAR>(), 10.0, 1e-12);

  EXPECT_NEAR(copy.get<parameters::TEST_SCALAR>(), 5.0, 1e-12);
  EXPECT_NEAR(copy.get<parameters::TEST_VECTOR>()[0], 5.0, 1e-12);
  EXPECT_NEAR(copy["reporter"].get<parameters::TEST_SCALAR>(), 50.0, 1e-12);

  // ... and the other direction
  parameters.set<parameters::TEST_SCALAR>(7.0);
  EXPECT_NEAR(copy.get<parameters::TEST_SCALAR>(), 5.0, 1e-12);
}

TEST(parametersD, copyAssignmentIsIndependent)
{
  OverrideGuard guard;
  Parameters parameters;
  Parameters copy;

  parameters.set<parameters::TEST_SCALAR>(1.0);
  copy.set<parameters::TEST_SCALAR>(2.0);

  copy = parameters;
  EXPECT_NEAR(copy.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);

  copy.set<parameters::TEST_SCALAR>(3.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
}

/// The accessor form evaluates against the set it is read through, so a copy
/// tracks the copy's own base rather than the original's
TEST(parametersD, accessorFormDependencyFollowsTheCopy)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(1.0);
  parameters.set<parameters::TEST_DERIVED>([](auto& p) {
    return 2.0 * p.template get<parameters::TEST_BASE>();
  });

  Parameters copy = parameters;
  copy.set<parameters::TEST_BASE>(10.0);

  EXPECT_NEAR(copy.get<parameters::TEST_DERIVED>(), 20.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 2.0, 1e-12);
}

TEST(parametersD, accessorFormTracksItsBase)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(1.0);
  parameters.set<parameters::TEST_DERIVED>([](auto& p) {
    return 2.0 * p.template get<parameters::TEST_BASE>();
  });
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 2.0, 1e-12);

  parameters.set<parameters::TEST_BASE>(4.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 8.0, 1e-12);
}

TEST(parametersD, accessorFormDetectsCycles)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>([](auto& p) {
    return p.template get<parameters::TEST_DERIVED>() + 1.0;
  });
  parameters.set<parameters::TEST_DERIVED>([](auto& p) {
    return p.template get<parameters::TEST_BASE>() + 1.0;
  });

  EXPECT_THROW(parameters.get<parameters::TEST_BASE>(), std::runtime_error);
}

/// Naming the type concretely avoids the dependent-name 'template' keyword
TEST(parametersD, accessorFormWithConcreteType)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_BASE>(1.0);
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 2.0 * p.get<parameters::TEST_BASE>();
  });

  Parameters copy = parameters;
  copy.set<parameters::TEST_BASE>(10.0);

  EXPECT_NEAR(copy.get<parameters::TEST_DERIVED>(), 20.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 2.0, 1e-12);
}

TEST(parametersD, moveKeepsValues)
{
  OverrideGuard guard;
  Parameters parameters;

  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters["reporter"].set<parameters::TEST_SCALAR>(10.0);

  Parameters moved = std::move(parameters);

  EXPECT_NEAR(moved.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
  EXPECT_NEAR(moved["reporter"].get<parameters::TEST_SCALAR>(), 10.0, 1e-12);
}

TEST(parametersD, moveAssignmentKeepsValues)
{
  OverrideGuard guard;
  Parameters parameters;
  Parameters target;

  parameters.set<parameters::TEST_SCALAR>(1.0);
  target = std::move(parameters);

  EXPECT_NEAR(target.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
}

TEST(parametersD, copyAsADf)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.5);
  parameters.set<parameters::TEST_VECTOR>(Vector<T,2>{1.0, 2.0});
  parameters["reporter"].set<parameters::TEST_SCALAR>(10.0);

  Parameters copy = parameters;

  auto adfParameters = copy.copyAs<ADf>();

  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 1.5, 1e-12);
  EXPECT_NEAR(adfParameters.get<parameters::TEST_VECTOR>()[1].v(), 2.0, 1e-12);
  EXPECT_NEAR(adfParameters["reporter"].get<parameters::TEST_SCALAR>().v(), 10.0, 1e-12);
}


TEST(parametersD, overlayConvertsOnAccess)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.5);
  parameters.set<parameters::TEST_VECTOR>(Vector<T,2>{1.0, 2.0});
  parameters.set<parameters::TEST_COUNT>(3);
  parameters["reporter"].set<parameters::TEST_OTHER>(10.0);

  auto adfParameters = parameters.copyAs<ADf>();

  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 1.5, 1e-12);
  EXPECT_NEAR(adfParameters.get<parameters::TEST_VECTOR>()[1].v(), 2.0, 1e-12);
  EXPECT_EQ(adfParameters.get<parameters::TEST_COUNT>(), 3);
  EXPECT_NEAR(adfParameters["reporter"].get<parameters::TEST_OTHER>().v(), 10.0, 1e-12);
}

/// The snapshot is owned, so the overlay outlives the instance it came from
TEST(parametersD, overlayOutlivesItsSource)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  ParametersD<ADf,DESCRIPTOR> adfParameters;
  {
    Parameters parameters;
    parameters.set<parameters::TEST_SCALAR>(2.5);
    adfParameters = parameters.copyAs<ADf>();
  }

  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 2.5, 1e-12);
}

/// Values set on the overlay shadow the snapshot
TEST(parametersD, overlayValuesShadowTheSource)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  auto adfParameters = parameters.copyAs<ADf>();
  adfParameters.set<parameters::TEST_SCALAR>(ADf(9.0));

  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 9.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
}

/// Later changes to the origin do not reach the overlay: it holds a snapshot
TEST(parametersD, overlayIsASnapshot)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  auto adfParameters = parameters.copyAs<ADf>();

  parameters.set<parameters::TEST_SCALAR>(7.0);

  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 1.0, 1e-12);
}

/// A chain double -> ADf<1> -> ADf<2> resolves through both links
TEST(parametersD, chainedOverlayResolves)
{
  OverrideGuard guard;
  using ADf1 = util::ADf<T,1>;
  using ADf2 = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(3.5);
  parameters["reporter"].set<parameters::TEST_OTHER>(4.5);

  auto first = parameters.copyAs<ADf1>();
  auto second = first.copyAs<ADf2>();

  // Nothing was ever read from first, so the value must come through both links
  EXPECT_NEAR(second.get<parameters::TEST_SCALAR>().v(), 3.5, 1e-12);
  EXPECT_NEAR(second["reporter"].get<parameters::TEST_OTHER>().v(), 4.5, 1e-12);
}

TEST(parametersD, chainedOverlayShadowsAtEachLink)
{
  OverrideGuard guard;
  using ADf1 = util::ADf<T,1>;
  using ADf2 = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  auto first = parameters.copyAs<ADf1>();
  first.set<parameters::TEST_SCALAR>(ADf1(5.0));

  auto second = first.copyAs<ADf2>();

  EXPECT_NEAR(second.get<parameters::TEST_SCALAR>().v(), 5.0, 1e-12);
}

/// Each link owns a copy, so a cycle of overlays terminates
TEST(parametersD, mutualOverlayTerminates)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  auto adfParameters = parameters.copyAs<ADf>();
  parameters = adfParameters.copyAs<T>();

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 1.0, 1e-12);
}

TEST(parametersD, overlayOfDependentParameterIsEvaluated)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_BASE>(2.0);
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 3.0 * p.get<parameters::TEST_BASE>();
  });

  auto adfParameters = parameters.copyAs<ADf>();

  EXPECT_NEAR(adfParameters.get<parameters::TEST_DERIVED>().v(), 6.0, 1e-12);
}

/// A copied overlay shares the snapshot rather than duplicating the chain
/// An overlay must not consult the instance it was made from after construction
TEST(parametersD, overlayOfDependentIgnoresLaterSourceChange)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_BASE>(2.0);
  parameters.set<parameters::TEST_DERIVED>([](auto& p) {
    return 3.0 * p.template get<parameters::TEST_BASE>();
  });

  auto adfParameters = parameters.copyAs<ADf>();

  // Changing the origin after the snapshot was taken must not reach the overlay
  parameters.set<parameters::TEST_BASE>(100.0);

  EXPECT_NEAR(adfParameters.get<parameters::TEST_DERIVED>().v(), 6.0, 1e-12);
}

/// Reading through an overlay must leave the instance it was made from alone
TEST(parametersD, overlayDoesNotMutateItsOrigin)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_BASE>(2.0);
  parameters.set<parameters::TEST_DERIVED>([](auto& p) {
    return 3.0 * p.template get<parameters::TEST_FILLER_0>()
         + 3.0 * p.template get<parameters::TEST_BASE>();
  });

  auto adfParameters = parameters.copyAs<ADf>();

  const std::size_t before = parameters.version();
  adfParameters.get<parameters::TEST_DERIVED>();
  EXPECT_EQ(parameters.version(), before);
}

/// A base set on the copy must feed the copy's derived parameters
TEST(parametersD, settingBaseOnCopyUpdatesDerived)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_BASE>(2.0);
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 3.0 * p.get<parameters::TEST_BASE>();
  });

  auto adfParameters = parameters.copyAs<ADf>();
  adfParameters.set<parameters::TEST_BASE>(ADf(10.0));

  EXPECT_NEAR(adfParameters.get<parameters::TEST_DERIVED>().v(), 30.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_DERIVED>(), 6.0, 1e-12);
}

/// Reads see the copy's own values, not the snapshot's
TEST(parametersD, copyReadsItsOwnValues)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  auto adfParameters = parameters.copyAs<ADf>();
  adfParameters.set<parameters::TEST_SCALAR>(ADf(5.0));

  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 5.0, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 1.0, 1e-12);
}

/// Setting a derived parameter on the copy stops it deriving there for good
TEST(parametersD, settingDerivedOnCopyOverridesEvaluation)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_BASE>(2.0);
  parameters.set<parameters::TEST_DERIVED>([](Parameters& p) {
    return 3.0 * p.get<parameters::TEST_BASE>();
  });

  auto adfParameters = parameters.copyAs<ADf>();
  adfParameters.set<parameters::TEST_DERIVED>(ADf(99.0));
  adfParameters.set<parameters::TEST_BASE>(ADf(10.0));

  EXPECT_NEAR(adfParameters.get<parameters::TEST_DERIVED>().v(), 99.0, 1e-12);
}

TEST(parametersD, copiedOverlayStillResolves)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.5);

  auto adfParameters = parameters.copyAs<ADf>();

  ParametersD<ADf,DESCRIPTOR> copy = adfParameters;
  EXPECT_NEAR(copy.get<parameters::TEST_SCALAR>().v(), 1.5, 1e-12);

  copy.set<parameters::TEST_SCALAR>(ADf(8.0));
  EXPECT_NEAR(adfParameters.get<parameters::TEST_SCALAR>().v(), 1.5, 1e-12);
}

/// An overlay must not invent values the source never had
TEST(parametersD, overlayFallsBackToDefaults)
{
  OverrideGuard guard;
  using ADf = util::ADf<T,2>;

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);

  auto adfParameters = parameters.copyAs<ADf>();

  EXPECT_FALSE(adfParameters.provides<parameters::TEST_OTHER>());
  EXPECT_NEAR(adfParameters.get<parameters::TEST_OTHER>().v(), 0.0, 1e-12);
}


TEST(parametersD, overrideFromCli)
{
  OverrideGuard guard;
  setOverrides({"--TEST_SCALAR", "3.5", "--TEST_COUNT", "7"});

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters.set<parameters::TEST_COUNT>(1);

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 3.5, 1e-12);
  EXPECT_EQ(parameters.get<parameters::TEST_COUNT>(), 7);
}

TEST(parametersD, overrideFromCliReachesSubgroups)
{
  OverrideGuard guard;
  setOverrides({"--reporter.vtk.TEST_SCALAR", "8.5"});

  Parameters parameters;
  parameters["reporter"]["vtk"].set<parameters::TEST_SCALAR>(1.0);

  EXPECT_NEAR(parameters["reporter"]["vtk"].get<parameters::TEST_SCALAR>(), 8.5, 1e-12);
}

TEST(parametersD, laterSetWinsOverCli)
{
  OverrideGuard guard;
  setOverrides({"--TEST_SCALAR", "3.5"});

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 3.5, 1e-12);

  parameters.set<parameters::TEST_SCALAR>(9.0);
  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 9.0, 1e-12);
}

TEST(parametersD, overrideFromXml)
{
  OverrideGuard guard;

  const std::string filename = "testParametersD_case.xml";
  {
    std::ofstream out(filename, std::ios::trunc);
    out << "<?xml version=\"1.0\"?>\n"
        << "<Case>\n<Parameters>\n"
        << "  <TEST_SCALAR>2.5</TEST_SCALAR>\n"
        << "  <reporter><vtk><TEST_SCALAR>4.5</TEST_SCALAR></vtk></reporter>\n"
        << "</Parameters>\n</Case>\n";
  }
  setOverrides({"--case-file", filename});

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters["reporter"]["vtk"].set<parameters::TEST_SCALAR>(1.0);

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 2.5, 1e-12);
  EXPECT_NEAR(parameters["reporter"]["vtk"].get<parameters::TEST_SCALAR>(), 4.5, 1e-12);

  std::remove(filename.c_str());
}

TEST(parametersD, cliWinsOverXml)
{
  OverrideGuard guard;

  const std::string filename = "testParametersD_case.xml";
  {
    std::ofstream out(filename, std::ios::trunc);
    out << "<?xml version=\"1.0\"?>\n"
        << "<Case>\n<Parameters>\n"
        << "  <TEST_SCALAR>2.5</TEST_SCALAR>\n"
        << "  <TEST_OTHER>6.5</TEST_OTHER>\n"
        << "</Parameters>\n</Case>\n";
  }
  setOverrides({"--case-file", filename, "--TEST_SCALAR", "3.5"});

  Parameters parameters;
  parameters.set<parameters::TEST_SCALAR>(1.0);
  parameters.set<parameters::TEST_OTHER>(1.0);

  EXPECT_NEAR(parameters.get<parameters::TEST_SCALAR>(), 3.5, 1e-12);
  EXPECT_NEAR(parameters.get<parameters::TEST_OTHER>(), 6.5, 1e-12);

  std::remove(filename.c_str());
}
