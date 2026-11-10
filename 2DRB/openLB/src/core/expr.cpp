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

#include <iostream>

#include "olbInit.h"
#include "expr.h"

#include "utilities/oalgorithm.h"

namespace olb {

std::array<std::size_t, static_cast<std::size_t>(Expr::Op::End)> Expr::_stats = { };

namespace {

thread_local std::size_t nextNodeId = 0;

struct BinaryKey {
  Expr::Op op;
  const void* lhs;
  const void* rhs;
  bool operator==(const BinaryKey& o) const {
    return op == o.op && lhs == o.lhs && rhs == o.rhs;
  }
};

struct UnaryKey {
  Expr::Op op;
  const void* arg;
  bool operator==(const UnaryKey& o) const {
    return op == o.op && arg == o.arg;
  }
};

struct BinaryKeyHash {
  std::size_t operator()(const BinaryKey& k) const {
    std::size_t h = std::hash<std::size_t>{}(static_cast<std::size_t>(k.op));
    h ^= std::hash<const void*>{}(k.lhs) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    h ^= std::hash<const void*>{}(k.rhs) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    return h;
  }
};

struct UnaryKeyHash {
  std::size_t operator()(const UnaryKey& k) const {
    std::size_t h = std::hash<std::size_t>{}(static_cast<std::size_t>(k.op));
    h ^= std::hash<const void*>{}(k.arg) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    return h;
  }
};

}

void Expr::describeNode(const Node& node, std::stringstream& out) {
  switch (node.kind) {
  case Node::Kind::Symbol:
    out << node.name;
    return;
  case Node::Kind::Constant:
    out << std::to_string(node.value);
    return;
  case Node::Kind::Binary:
    switch (node.op) {
    case Op::Add: out << "("; describeNode(*node.lhs, out); out << ") + ("; describeNode(*node.rhs, out); out << ")"; return;
    case Op::Sub: out << "("; describeNode(*node.lhs, out); out << ") - ("; describeNode(*node.rhs, out); out << ")"; return;
    case Op::Mul: out << "("; describeNode(*node.lhs, out); out << ") * ("; describeNode(*node.rhs, out); out << ")"; return;
    case Op::Div: out << "("; describeNode(*node.lhs, out); out << ") / ("; describeNode(*node.rhs, out); out << ")"; return;
    case Op::Pow: out << "Pow("; describeNode(*node.lhs, out); out << ","; describeNode(*node.rhs, out); out << ")"; return;
    default: throw std::invalid_argument("Unsupported binary operation");
    }
  case Node::Kind::Unary:
    switch (node.op) {
    case Op::Sqrt: out << "sqrt("; describeNode(*node.lhs, out); out << ")"; return;
    case Op::Abs:  out << "Abs(";  describeNode(*node.lhs, out); out << ")"; return;
    case Op::Exp:  out << "exp(";  describeNode(*node.lhs, out); out << ")"; return;
    default: throw std::invalid_argument("Unsupported unary operation");
    }
  }
}

std::size_t Expr::nodeSize(const Node& node) {
  switch (node.kind) {
  case Node::Kind::Symbol:
  case Node::Kind::Constant:
    return 1;
  case Node::Kind::Unary:
    return 1 + nodeSize(*node.lhs);
  case Node::Kind::Binary:
    return 1 + nodeSize(*node.lhs) + nodeSize(*node.rhs);
  }
  return 1;
}

Expr::Ref Expr::internConstant(double value) {
  static thread_local std::unordered_map<double, Ref> table;
  auto [it, inserted] = table.try_emplace(value, nullptr);
  if (inserted) {
    auto node = std::make_shared<Node>();
    node->kind = Node::Kind::Constant;
    node->value = value;
    node->id = nextNodeId++;
    it->second = node;
  }
  return it->second;
}

Expr::Ref Expr::internSymbol(const std::string& name) {
  static thread_local std::unordered_map<std::string, Ref> table;
  auto [it, inserted] = table.try_emplace(name, nullptr);
  if (inserted) {
    auto node = std::make_shared<Node>();
    node->kind = Node::Kind::Symbol;
    node->name = name;
    node->id = nextNodeId++;
    it->second = node;
  }
  return it->second;
}

Expr::Ref Expr::internBinary(Op op, Ref lhs, Ref rhs) {
  static thread_local std::unordered_map<BinaryKey, Ref, BinaryKeyHash> table;
  // Order commutative operands so a+b and b+a share a node (swap is bit-exact;
  // associativity is never altered).
  if (op == Op::Add || op == Op::Mul) {
    if (rhs->id < lhs->id) {
      std::swap(lhs, rhs);
    }
  }
  BinaryKey key{op, lhs.get(), rhs.get()};
  auto [it, inserted] = table.try_emplace(key, nullptr);
  if (inserted) {
    auto node = std::make_shared<Node>();
    node->kind = Node::Kind::Binary;
    node->op = op;
    node->lhs = std::move(lhs);
    node->rhs = std::move(rhs);
    node->id = nextNodeId++;
    it->second = node;
  }
  return it->second;
}

Expr::Ref Expr::internUnary(Op op, Ref arg) {
  static thread_local std::unordered_map<UnaryKey, Ref, UnaryKeyHash> table;
  UnaryKey key{op, arg.get()};
  auto [it, inserted] = table.try_emplace(key, nullptr);
  if (inserted) {
    auto node = std::make_shared<Node>();
    node->kind = Node::Kind::Unary;
    node->op = op;
    node->lhs = std::move(arg);
    node->id = nextNodeId++;
    it->second = node;
  }
  return it->second;
}

Expr::Expr(Ref node):
  _node(std::move(node)) { }

Expr::Expr():
  _node(internConstant(0.0)) { }

Expr::Expr(Expr&& rhs):
  _node(std::move(rhs._node)) { }

Expr::Expr(const Expr& rhs):
  _node(rhs._node) { }

Expr::Expr(double v):
  _node(internConstant(v)) { }

Expr::Expr(std::string name):
  _node(internSymbol(name)) { }

Expr::Expr(Expr lhs, Op op, Expr rhs):
  _node(internBinary(op, lhs._node, rhs._node)) { }

Expr::Expr(Op op, Expr rhs):
  _node(internUnary(op, rhs._node)) { }

Expr& Expr::operator=(const Expr& rhs) {
  _node = rhs._node;
  return *this;
}

Expr& Expr::operator+=(Expr rhs) {
  if (rhs.isZero()) return *this;
  if (isZero()) { _node = rhs._node; return *this; }

  increment(Expr::Op::Add);
  _node = internBinary(Op::Add, _node, rhs._node);
  return *this;
}

Expr& Expr::operator-=(Expr rhs) {
  if (rhs.isZero()) return *this;

  increment(Expr::Op::Sub);
  _node = internBinary(Op::Sub, _node, rhs._node);
  return *this;
}

Expr& Expr::operator*=(Expr rhs) {
  if (rhs.isZero() || isZero()) {
    _node = internConstant(0.0);
    return *this;
  }

  if (rhs.isOne()) return *this;

  increment(Expr::Op::Mul);
  _node = internBinary(Op::Mul, _node, rhs._node);
  return *this;
}

Expr& Expr::operator/=(Expr rhs) {
  increment(Expr::Op::Div);
  _node = internBinary(Op::Div, _node, rhs._node);
  return *this;
}

Expr::operator int() const {
  throw std::domain_error("Expr type can not be converted to int");
}

Expr::operator long unsigned int() const {
  throw std::domain_error("Expr type can not be converted to long unsigned int");
}

void Expr::describe(std::stringstream& out) const {
  describeNode(*_node, out);
};

std::string Expr::describe() const {
  std::stringstream out;
  describe(out);
  return out.str();
};

bool Expr::isSymbol(std::string name) const {
  return _node->kind == Node::Kind::Symbol && _node->name == name;
}

std::size_t Expr::size() const {
  return nodeSize(*_node);
}

Expr operator+(Expr lhs, Expr rhs) {
  if (lhs.isZero()) return rhs;   // 0 + x = x
  if (rhs.isZero()) return lhs;   // x + 0 = x

  Expr::increment(Expr::Op::Add);
  return Expr(lhs, Expr::Op::Add, rhs);
}

Expr operator-(Expr lhs, Expr rhs) {
  if (rhs.isZero()) return lhs;   // x - 0 = x

  Expr::increment(Expr::Op::Sub);
  return Expr(lhs, Expr::Op::Sub, rhs);
}

Expr operator*(Expr lhs, Expr rhs) {
  if (lhs.isZero() || rhs.isZero()) return Expr(0.0); // x*0 or 0*x = 0
  if (lhs.isOne()) return rhs; // 1*x = x
  if (rhs.isOne()) return lhs; // x*1 = x

  Expr::increment(Expr::Op::Mul);
  return Expr(lhs, Expr::Op::Mul, rhs);
}

Expr operator/(Expr lhs, Expr rhs) {
  if (lhs.isZero()) return Expr(0.0); // 0/x = 0
  Expr::increment(Expr::Op::Div);
  return Expr(lhs, Expr::Op::Div, rhs);
}

Expr operator-(Expr rhs) {
  return Expr(-1.) * rhs;
}

Expr operator%(Expr lhs, int rhs) {
  throw std::domain_error("Invalid operator for Expr type: '%'");
}

bool operator==(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '=='");
}

bool operator!=(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '!='");
}

bool operator>(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '>'");
}

bool operator<(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '<'");
}

bool operator>=(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '>='");
}

bool operator<=(const Expr& rhs, const Expr& lhs) {
  throw std::domain_error("Invalid operator for Expr type: '<='");
}


void Expr::DagContext::countEdges(const Node* node) {
  // Visit each unique node once, count one in-edge per child reference.
  if (!_seen.insert(node).second) {
    return;
  }
  if (node->kind == Node::Kind::Binary) {
    _indegree[node->lhs.get()]++;
    _indegree[node->rhs.get()]++;
    countEdges(node->lhs.get());
    countEdges(node->rhs.get());
  } else if (node->kind == Node::Kind::Unary) {
    _indegree[node->lhs.get()]++;
    countEdges(node->lhs.get());
  }
}

void Expr::DagContext::scan(const Expr& expr) {
  _indegree[expr._node.get()]++; // count the root reference
  countEdges(expr._node.get());
}

std::string Expr::DagContext::buildRhs(const Node* node,
                                       const std::string& lhs, const std::string& rhs) {
  switch (node->kind) {
  case Node::Kind::Binary:
    switch (node->op) {
    case Op::Add: return "(" + lhs + ") + (" + rhs + ")";
    case Op::Sub: return "(" + lhs + ") - (" + rhs + ")";
    case Op::Mul: return "(" + lhs + ") * (" + rhs + ")";
    case Op::Div: return "(" + lhs + ") / (" + rhs + ")";
    case Op::Pow: return "Pow(" + lhs + "," + rhs + ")";
    default: throw std::invalid_argument("Unsupported binary operation");
    }
  case Node::Kind::Unary:
    switch (node->op) {
    case Op::Sqrt: return "sqrt(" + lhs + ")";
    case Op::Abs:  return "Abs(" + lhs + ")";
    case Op::Exp:  return "exp(" + lhs + ")";
    default: throw std::invalid_argument("Unsupported unary operation");
    }
  default:
    throw std::invalid_argument("Leaf node has no composite right-hand side");
  }
}

std::string Expr::DagContext::emit(const Node* node) {
  if (node->kind == Node::Kind::Symbol) {
    return node->name;
  }
  if (node->kind == Node::Kind::Constant) {
    return std::to_string(node->value);
  }

  const bool shared = _indegree[node] >= 2;
  if (shared) {
    auto it = _emitted.find(node);
    if (it != _emitted.end()) {
      return it->second;
    }
  }

  std::string lhs = emit(node->lhs.get());
  std::string rhs = (node->kind == Node::Kind::Binary) ? emit(node->rhs.get()) : std::string{};
  std::string body = buildRhs(node, lhs, rhs);

  if (shared) {
    std::string sym = "t" + std::to_string(_counter++);
    std::string symref = "Symbol(\"" + sym + "\")";
    _bindings.push_back("Assignment(" + symref + ", " + body + ")");
    _emitted[node] = symref;
    return symref;
  }
  return body;
}

std::string Expr::DagContext::ref(const Expr& expr) {
  return emit(expr._node.get());
}

namespace util {

Expr sqrt(Expr x) {
  Expr::increment(Expr::Op::Sqrt);
  return Expr(Expr::Op::Sqrt, x);
}

Expr fabs(Expr x) {
  Expr::increment(Expr::Op::Abs);
  return Expr(Expr::Op::Abs, x);
}

Expr pow(Expr base, Expr exp) {
  Expr::increment(Expr::Op::Pow);
  return Expr(base, Expr::Op::Pow, exp);
}

Expr exp(Expr x) {
  Expr::increment(Expr::Op::Exp);
  return Expr(Expr::Op::Exp, x);
}

Expr max(Expr a, Expr b) {
  throw std::domain_error("Invalid operator for Expr type: 'max'");
}

Expr min(Expr a, Expr b) {
  throw std::domain_error("Invalid operator for Expr type: 'min'");
}

}

}
