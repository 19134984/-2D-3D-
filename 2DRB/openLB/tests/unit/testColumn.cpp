/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Clara Schragmann
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

using namespace olb;
typedef double T;


TEST(Column, basicTest)
{
  Column<T> c(100);
  for (std::size_t i = 0; i < 100; ++i) {
    c[i] = i;
  }
  for (std::size_t i = 0; i < 100; ++i) {
    EXPECT_EQ(c[i], i);
  }
}

TEST(Column, testSerializableSize)
{
  Column<T> a(100);

  for (std::size_t i = 0; i < 100; ++i) {
    a[i] = i;
  }

  EXPECT_EQ((100*sizeof(T)+sizeof(std::size_t)),a.getSerializableSize());
}

TEST(Column, testData)
{
  Column<T> a(100);
  Column<T> b(100);

  for (std::size_t i = 0; i < 100; ++i) {
    a[i] = i;
  }

  a.save("test.dat");
  b.load("test.dat");

  for (std::size_t i = 0; i < 100; ++i) {
    EXPECT_EQ(a[i], b[i]);
  }
}

TEST(Column, resize)
{
  Column<T> c(100);
  for (std::size_t i = 0; i < 100; ++i) {
    c[i] = i;
  }
  c.resize(200);
  for (std::size_t i = 0; i < 100; ++i) {
    EXPECT_EQ(c[i], i);
  }
  c.resize(50);
  for (std::size_t i = 0; i < 50; ++i) {
    EXPECT_EQ(c[i], i);
  }
}

TEST(ColumnVector, basic)
{
  ColumnVector<Column<T>,5> a(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
    }
  }

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      EXPECT_EQ((a[i][j]),(i+j));
    }
  }
}

TEST(ColumnVector, VectorPtr)
{
  ColumnVector<Column<T>,5> a(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row = a.getRowPointer(j);
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_EQ(row[i], i+j);
    }
  }
}

TEST(ColumnVector, VectorPtrAdd)
{
  ColumnVector<Column<T>,5> a(100);
  ColumnVector<Column<T>,5> b(100);
  Vector<T,5> c{1,2,3,4,5};

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
      (b[i])[j] = i + j;
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row = a.getRowPointer(j);
    row+=2.;
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_NEAR(row[i], i+j+2., 1e-10);
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    row_a+=row_b;
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_NEAR(row_a[i], i+j+i+j+2., 1e-10);
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_b2 = b.getRowPointer(j);
    row_b2+=c;
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_NEAR(row_b2[i], i+j+i+1., 1e-10);
    }
  }
}

TEST(ColumnVector, VectorPtrSubtract)
{
  ColumnVector<Column<T>,5> a(100);
  ColumnVector<Column<T>,5> b(100);
  Vector<T,5> c{1,2,3,4,5};

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = j + i;
      (b[i])[j] =     i;
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row = a.getRowPointer(j);
    row-=2.;
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_NEAR(row[i], i+j-2., 1e-10);
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    row_a-=row_b;
    for (std::size_t i=0; i < 5; ++i) {
      EXPECT_NEAR(row_a[i], i+j-2.-i, 1e-10);
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_b2 = b.getRowPointer(j);
    row_b2-=c;
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_NEAR(row_b2[i], i-i-1., 1e-10);
    }
  }
}

TEST(ColumnVector, VectorPtrMultiplyAndDivide)
{
  ColumnVector<Column<T>,5> a(100);
  ColumnVector<Column<T>,5> b(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] =    i + j;
      (b[i])[j] = 2*(i + j);
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    double scalarproduct=0;
    for (std::size_t i = 0; i < 5; i++) {
      scalarproduct += row_a[i]*row_b[i];
    }
    EXPECT_NEAR(row_a*row_b, scalarproduct, 1e-10);
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    row_a*=2.;
    row_b/=2.;
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_NEAR(row_a[i], (i+j)*2, 1e-10);
      EXPECT_NEAR(row_b[i], (i+j), 1e-10);
    }
  }
}

TEST(ColumnVector, VectorPtrComparisonEquals)
{
  ColumnVector<Column<T>,5> a(100);
  ColumnVector<Column<T>,5> b(100);
  ColumnVector<Column<T>,5> c(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
      (b[i])[j] = i + j;
      (c[i])[j] = i + j + 5.;
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    auto row_c = c.getRowPointer(j);
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_TRUE(row_a[i]==row_b[i]);
      EXPECT_FALSE(row_a[i]==row_c[i]);
      EXPECT_TRUE(row_a[i]!=row_c[i]);
      EXPECT_FALSE(row_a[i]!=row_b[i]);
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    auto row_c = c.getRowPointer(j);

    EXPECT_TRUE(row_a==row_b);
    EXPECT_FALSE(row_a==row_c);
    EXPECT_TRUE(row_a!=row_c);
    EXPECT_FALSE(row_a!=row_b);
  }
}

TEST(ColumnVector, VectorPtrComparison)
{
  ColumnVector<Column<T>,5> a(100);
  ColumnVector<Column<T>,5> b(100);
  ColumnVector<Column<T>,5> c(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
      (b[i])[j] = i + j;
      (c[i])[j] = i + j + 5.;
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    auto row_b = b.getRowPointer(j);
    auto row_c = c.getRowPointer(j);

    EXPECT_TRUE(row_a<=row_b);
    EXPECT_FALSE(row_a<row_b);
    EXPECT_TRUE(row_b>=row_a);
    EXPECT_FALSE(row_b>row_a);
    EXPECT_TRUE(row_b<=row_c);
    EXPECT_TRUE(row_b<row_c);
    EXPECT_TRUE(row_c>=row_b);
    EXPECT_TRUE(row_c>row_b);
  }
}

TEST(ColumnVector, testValueReassigning)
{
  ColumnVector<Column<T>,5> a(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
    }
  }

  for (std::size_t j = 0; j < 100; ++j) {
    auto row_a = a.getRowPointer(j);
    for (std::size_t i = 0; i < 5; ++i) {
      row_a[i]= i + j + 5;
    }
  }

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      EXPECT_EQ((a[i])[j],i+j+5);
    }
  }
}

TEST(ColumnVector, testSerializableSize)
{
  ColumnVector<Column<T>,5> a(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
    }
  }

  EXPECT_EQ((5*(100*sizeof(T)+sizeof(std::size_t))),a.getSerializableSize());
}

TEST(ColumnVector, testData)
{
  ColumnVector<Column<T>,5> a(100);
  ColumnVector<Column<T>,5> b(100);

  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 100; ++j) {
      (a[i])[j] = i + j;
    }
  }
  a.save("test.dat");
  b.load("test.dat");

  for (std::size_t j = 0; j < 100; ++j) {
    for (std::size_t i = 0; i < 5; ++i) {
      EXPECT_EQ(a[i][j],b[i][j]);
    }
  }
}
