/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CORE_EXPR_H
#define CORE_EXPR_H

#include <string>
#include <sstream>
#include <memory>
#include <array>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <numeric>

#include "core/meta.h"

namespace olb {

/// Extraction of expression trees for code generation
/**
 * Expr is hash-consed: structurally identical subexpressions are interned into a
 * single shared, immutable node, so the in-memory representation is a DAG.
 * Addition and multiplication canonically order their operands, so commutatively
 * equal subexpressions (a+b == b+a) share the same node. This means CSE is largely
 * performed during construction.
 */
class Expr : public ExprBase {
public:
  enum struct Op : std::size_t {
    Add, Mul, Sub, Div, Sqrt, Abs, Pow, Exp,
    End // keep last
  };

private:
  static std::array<std::size_t, static_cast<std::size_t>(Op::End)> _stats;

  /// Immutable, interned expression node; children point at other interned nodes.
  struct Node {
    enum struct Kind : std::uint8_t { Constant, Symbol, Binary, Unary };

    Kind kind;
    Op op = Op::End;                  // Binary / Unary
    double value = 0.0;               // Constant
    std::string name;                 // Symbol
    std::shared_ptr<const Node> lhs;  // Binary lhs / Unary argument
    std::shared_ptr<const Node> rhs;  // Binary rhs
    std::size_t id = 0;               // monotonic creation id (deterministic ordering)
  };

  using Ref = std::shared_ptr<const Node>;

  Ref _node;

  explicit Expr(Ref node);

  static Ref internConstant(double value);
  static Ref internSymbol(const std::string& name);
  static Ref internBinary(Op op, Ref lhs, Ref rhs);
  static Ref internUnary(Op op, Ref arg);

  static void describeNode(const Node& node, std::stringstream& out);
  static std::size_t nodeSize(const Node& node);

public:
  Expr();
  Expr(Expr&& rhs);
  Expr(const Expr& rhs);
  Expr(double v);
  Expr(std::string name);
  Expr(Expr lhs, Op op, Expr rhs);
  Expr(Op op, Expr rhs);

  Expr& operator=(const Expr& rhs);

  Expr& operator+=(Expr rhs);
  Expr& operator-=(Expr rhs);
  Expr& operator*=(Expr rhs);
  Expr& operator/=(Expr rhs);

  explicit operator int() const;
  explicit operator long unsigned int() const;

  /// Writes the serialized expression tree to out
  void describe(std::stringstream& out) const;
  /// Returns the serialized expression tree
  std::string describe() const;

  /// Helper function to quickly check if Expr is a specific symbol
  bool isSymbol(std::string name) const;

  /// Returns expanded size of tree
  std::size_t size() const;

  /// Resets FLOP counter
  static void reset() {
    _stats.fill(0);
  }

  /// Increment FLOP counter
  static void increment(Op op) {
    ++_stats[static_cast<std::size_t>(op)];
  }

  /// Return accumulated FLOPs
  static std::size_t count() {
    return std::accumulate(_stats.cbegin(), _stats.cend(), 0);
  }

  /// Return accumulated FLOPs of type
  static std::size_t count(Op op) {
    return _stats[static_cast<std::size_t>(op)];
  }

  bool isZero() const {
    return _node->kind == Node::Kind::Constant && _node->value == 0.0;
  }

  bool isOne() const {
    return _node->kind == Node::Kind::Constant && _node->value == 1.0;
  }

  /// Shared-subexpression serializer for the hash-consed DAG.
  /**
   * Usage: scan() every root expression, then ref() each root to obtain its
   * serialized right-hand side. Subexpressions referenced more than once are
   * emitted exactly once as `t<k>` bindings (retrievable via bindings()) and
   * referenced by symbol. Single-use subexpressions are inlined.
   */
  class DagContext {
  private:
    std::unordered_set<const Node*> _seen;
    std::unordered_map<const Node*, std::size_t> _indegree;
    std::unordered_map<const Node*, std::string> _emitted;
    std::vector<std::string> _bindings;
    std::size_t _counter = 0;

    void countEdges(const Node* node);
    static std::string buildRhs(const Node* node,
                                const std::string& lhs, const std::string& rhs);
    std::string emit(const Node* node);

  public:
    /// Record a root so shared subexpressions below it can be detected
    void scan(const Expr& expr);
    /// Serialize expr, emitting shared subexpressions as bindings; returns its rhs
    std::string ref(const Expr& expr);
    /// Bindings (`Assignment(...)` strings) collected so far, in dependency order
    const std::vector<std::string>& bindings() const {
      return _bindings;
    }
  };

};

Expr operator+(Expr lhs, Expr rhs);
Expr operator-(Expr lhs, Expr rhs);
Expr operator*(Expr lhs, Expr rhs);
Expr operator/(Expr lhs, Expr rhs);
Expr operator-(Expr rhs);

Expr operator%(Expr lhs, int rhs);

bool operator==(const Expr& rhs, const Expr& lhs);
bool operator!=(const Expr& rhs, const Expr& lhs);
bool operator>(const Expr& rhs, const Expr& lhs);
bool operator<(const Expr& rhs, const Expr& lhs);
bool operator>=(const Expr& rhs, const Expr& lhs);
bool operator<=(const Expr& rhs, const Expr& lhs);

namespace util {

Expr sqrt(Expr x);
Expr fabs(Expr x);
Expr pow(Expr base, Expr exp);
Expr exp(Expr x);


}

}

#endif
