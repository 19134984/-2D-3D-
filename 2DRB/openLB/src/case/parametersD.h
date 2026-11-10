/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender, Shota Ito
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

#ifndef CASE_PARAMETERSD_H
#define CASE_PARAMETERSD_H

#include "io/input.h"
#include "io/textTable.h"
#include "parameters.h"
#include "core/fieldParametersD.h"
#include "core/meta.h"

#include <functional>
#include <stdexcept>
#include <memory>
#include <optional>
#include <sstream>
#include <algorithm>
#include <map>
#include <vector>
#include <array>

namespace olb {

/// On-demand allocating parameter field storage
/// Precedence set<>() < XML < CLI, resolved on first get<>(); a later set<>() wins
template <concepts::BaseType T, concepts::Descriptor DESCRIPTOR>
class ParametersD final {
private:
  template <concepts::BaseType, concepts::Descriptor> friend class ParametersD;

  /// Cloneable owner of a single parameter value
  struct AbstractFieldD {
    virtual ~AbstractFieldD() = default;
    virtual std::unique_ptr<AbstractFieldD> clone() const = 0;
  };

  template <concepts::Field FIELD>
  struct ConcreteFieldD final : public AbstractFieldD {
    FieldD<T,DESCRIPTOR,FIELD> data;

    ConcreteFieldD(const FieldD<T,DESCRIPTOR,FIELD>& value)
      : data{value}
    { }

    std::unique_ptr<AbstractFieldD> clone() const override {
      return std::make_unique<ConcreteFieldD>(*this);
    }
  };

  /// Type-erased owner of a snapshot of parameters of a different base type
  struct AbstractSource {
    virtual ~AbstractSource() = default;
    virtual std::unique_ptr<AbstractSource> clone() const = 0;
  };

  template <concepts::BaseType S>
  struct ConcreteSource final : public AbstractSource {
    ParametersD<S,DESCRIPTOR> data;

    ConcreteSource(const ParametersD<S,DESCRIPTOR>& source)
      : data{source}
    { }

    std::unique_ptr<AbstractSource> clone() const override {
      return std::make_unique<ConcreteSource>(*this);
    }
  };

  /// Reads the components of one field out of a snapshot, in the requested scalar type
  using FetchSourceF = meta::function<bool(AbstractSource&,
                                          const std::vector<std::string>& path,
                                          std::type_index,
                                          void* out,
                                          const std::type_info& requested,
                                          std::size_t n)>;

  struct TypeErasedFieldD {
    /// Absent until first evaluation
    std::unique_ptr<AbstractFieldD> field;

    std::string name;
    bool isNative = false;
    bool isFlag = false;

    meta::function<std::string(ParametersD&)> value;
    meta::function<void(ParametersD&, const std::string&)> parse;
    meta::function<bool(ParametersD&, void*, const std::type_info&, std::size_t)> exportComponents;
    meta::function<void(ParametersD&, const void*, const std::type_info&, std::size_t)> importComponents;
    meta::function<void(ParametersD&)> calculate;

    bool isBeingComputedRightNow = false;
    bool resolved = false;
    /// Cached from the snapshot, hence droppable
    bool fromSource = false;

    TypeErasedFieldD() = default;
    TypeErasedFieldD(TypeErasedFieldD&&) = default;
    TypeErasedFieldD& operator=(TypeErasedFieldD&&) = default;

    TypeErasedFieldD(const TypeErasedFieldD& rhs)
      : field{rhs.field ? rhs.field->clone() : nullptr}
      , name{rhs.name}
      , isNative{rhs.isNative}
      , isFlag{rhs.isFlag}
      , value{rhs.value}
      , parse{rhs.parse}
      , exportComponents{rhs.exportComponents}
      , importComponents{rhs.importComponents}
      , calculate{rhs.calculate}
      , resolved{rhs.resolved}
      , fromSource{rhs.fromSource}
    { }

    TypeErasedFieldD& operator=(const TypeErasedFieldD& rhs) {
      TypeErasedFieldD copy{rhs};
      *this = std::move(copy);
      return *this;
    }
  };

  std::unordered_map<std::type_index, TypeErasedFieldD> _map;

  /// Dot-path of this (sub)group ("" for root, "fluid", "fluid.wall")
  std::string _prefix;
  /// Enclosing group, null at the root
  ParametersD* _parent = nullptr;
  std::map<std::string, std::unique_ptr<ParametersD>> _children;
  std::size_t _version = 0;

  /// Writes one field back into a snapshot, so later evaluations there see it
  using PushSourceF = meta::function<void(AbstractSource&,
                                         const std::vector<std::string>& path,
                                         std::type_index,
                                         const void* in,
                                         const std::type_info& provided,
                                         std::size_t n)>;

  /// Snapshot of a differently typed parameter set, held by the outermost group only
  std::unique_ptr<AbstractSource> _source;
  FetchSourceF _fetch;
  PushSourceF _push;
  /// Guards against pushing a value back that was just read out of the snapshot
  bool _materializing = false;

  static std::vector<std::string> splitPrefix(const std::string& prefix) {
    std::vector<std::string> segments;
    std::size_t start = 0;
    while (start <= prefix.size() && !prefix.empty()) {
      const std::size_t dot = prefix.find('.', start);
      if (dot == std::string::npos) {
        segments.push_back(prefix.substr(start));
        break;
      }
      segments.push_back(prefix.substr(start, dot - start));
      start = dot + 1;
    }
    return segments;
  }

  template <concepts::Field FIELD>
  void setupTypeErasedField(TypeErasedFieldD& typeErasedField) {
    typeErasedField.name = parameters::name<FIELD>();
    typeErasedField.isNative = parameters::isNativeName(typeErasedField.name);
    typeErasedField.isFlag =
      std::is_same_v<typename FIELD::template value_type<T>, bool>;
    typeErasedField.value = [](ParametersD& self) -> std::string {
      std::stringstream out;
      out << *self.template getFieldPtr<FIELD>();
      return out.str();
    };
    typeErasedField.parse = [](ParametersD& self, const std::string& str) {
      if (auto v = FieldD<T,DESCRIPTOR,FIELD>::fromString(str)) {
        self.template set<FIELD>(*v);
      } else {
        throw std::invalid_argument("Cannot parse value '" + str + "' for parameter '"
                                    + std::string(parameters::name<FIELD>()) + "'");
      }
    };
    typeErasedField.exportComponents = [](ParametersD& self, void* out,
                                          const std::type_info& requested, std::size_t n) -> bool {
      using VALUE_TYPE = typename FIELD::template value_type<T>;
      constexpr std::size_t size = DESCRIPTOR::template size<FIELD>();
      if (n != size) {
        return false;
      }
      auto& value = *self.template getFieldPtr<FIELD>();
      // Values travel as plain scalars, so an ADf source contributes values, not derivatives
      if constexpr (std::is_floating_point_v<VALUE_TYPE> || util::is_adf_v<VALUE_TYPE>) {
        if (requested == typeid(BaseType<T>)) {
          auto* buffer = static_cast<BaseType<T>*>(out);
          for (std::size_t i=0; i < size; ++i) {
            buffer[i] = value[i];
          }
          return true;
        }
      }
      if (requested == typeid(VALUE_TYPE)) {
        auto* buffer = static_cast<VALUE_TYPE*>(out);
        for (std::size_t i=0; i < size; ++i) {
          buffer[i] = value[i];
        }
        return true;
      }
      return false;
    };
    typeErasedField.importComponents = [](ParametersD& self, const void* in,
                                          const std::type_info& provided, std::size_t n) {
      using VALUE_TYPE = typename FIELD::template value_type<T>;
      constexpr std::size_t size = DESCRIPTOR::template size<FIELD>();
      if (n != size) {
        return;
      }
      FieldD<T,DESCRIPTOR,FIELD> value;
      if constexpr (std::is_floating_point_v<VALUE_TYPE> || util::is_adf_v<VALUE_TYPE>) {
        if (provided == typeid(BaseType<T>)) {
          auto* buffer = static_cast<const BaseType<T>*>(in);
          for (std::size_t i=0; i < size; ++i) {
            value[i] = VALUE_TYPE(buffer[i]);
          }
          self.template set<FIELD>(value);
          return;
        }
      }
      if (provided == typeid(VALUE_TYPE)) {
        auto* buffer = static_cast<const VALUE_TYPE*>(in);
        for (std::size_t i=0; i < size; ++i) {
          value[i] = buffer[i];
        }
        self.template set<FIELD>(value);
      }
    };
  }

  /// Accept one field into the (sub)group at path, creating the group if needed
  void importAt(const std::vector<std::string>& path, std::type_index field,
                const void* in, const std::type_info& provided, std::size_t n) {
    ParametersD* node = this;
    for (const auto& segment : path) {
      node = &(*node)[segment];
    }
    if (auto entry = node->_map.find(field);
        entry != node->_map.end() && entry->second.importComponents) {
      entry->second.importComponents(*node, in, provided, n);
    }
    else if (_source && _push) {
      _push(*_source, path, field, in, provided, n);
    }
  }

  /// Serve one field of the (sub)group at path, chaining into our own source if needed
  bool exportAt(const std::vector<std::string>& path, std::type_index field,
                void* out, const std::type_info& requested, std::size_t n) {
    if (ParametersD* node = find(path)) {
      if (auto entry = node->_map.find(field);
          entry != node->_map.end() && entry->second.exportComponents) {
        return entry->second.exportComponents(*node, out, requested, n);
      }
    }
    // Not held here, but this snapshot may resolve from a further one
    if (_source && _fetch) {
      return _fetch(*_source, path, field, out, requested, n);
    }
    return false;
  }

  /// Pull FIELD out of the source snapshot, converting into our own value type
  template <concepts::Field FIELD>
  bool tryMaterializeFromSource() {
    ParametersD& owner = outermost();
    if (!owner._source || !owner._fetch) {
      return false;
    }
    using VALUE_TYPE = typename FIELD::template value_type<T>;
    constexpr std::size_t size = DESCRIPTOR::template size<FIELD>();
    const auto path = splitPrefix(_prefix);
    FieldD<T,DESCRIPTOR,FIELD> materialized;
    if constexpr (std::is_floating_point_v<VALUE_TYPE> || util::is_adf_v<VALUE_TYPE>) {
      std::array<BaseType<T>,size> buffer;
      if (!owner._fetch(*owner._source, path, typeid(FIELD), buffer.data(), typeid(BaseType<T>), size)) {
        return false;
      }
      for (std::size_t i=0; i < size; ++i) {
        materialized[i] = VALUE_TYPE(buffer[i]);
      }
    } else {
      if (!owner._fetch(*owner._source, path, typeid(FIELD), &materialized, typeid(VALUE_TYPE), size)) {
        return false;
      }
    }
    _materializing = true;
    try {
      set<FIELD>(materialized);
    } catch (...) {
      _materializing = false;
      throw;
    }
    _materializing = false;
    _map[typeid(FIELD)].fromSource = true;
    return true;
  }

  /// Hand a locally set value to the snapshot so fresh evaluations there agree
  template <concepts::Field FIELD>
  void pushToSource(ParametersD& owner, const FieldD<T,DESCRIPTOR,FIELD>& value) {
    using VALUE_TYPE = typename FIELD::template value_type<T>;
    constexpr std::size_t size = DESCRIPTOR::template size<FIELD>();
    const auto path = splitPrefix(_prefix);
    if constexpr (std::is_floating_point_v<VALUE_TYPE> || util::is_adf_v<VALUE_TYPE>) {
      std::array<BaseType<T>,size> buffer;
      for (std::size_t i=0; i < size; ++i) {
        buffer[i] = value[i];
      }
      owner._push(*owner._source, path, typeid(FIELD), buffer.data(), typeid(BaseType<T>), size);
    } else {
      owner._push(*owner._source, path, typeid(FIELD), &value, typeid(VALUE_TYPE), size);
    }
  }

  /// Drop values cached from the snapshot so they are fetched afresh
  void dropSourceValues() {
    for (auto& [typeId, typeErasedField] : _map) {
      if (typeErasedField.fromSource) {
        typeErasedField.field.reset();
        typeErasedField.fromSource = false;
        typeErasedField.resolved = false;
      }
    }
    for (auto& [name, child] : _children) {
      child->dropSourceValues();
    }
  }

  /// The (sub)group at path, or null when it was never created
  ParametersD* find(const std::vector<std::string>& path) {
    ParametersD* node = this;
    for (const auto& segment : path) {
      auto iter = node->_children.find(segment);
      if (iter == node->_children.end()) {
        return nullptr;
      }
      node = iter->second.get();
    }
    return node;
  }

  void adoptChildren() {
    for (auto& [name, child] : _children) {
      child->_parent = this;
    }
  }

  void attachSource(std::unique_ptr<AbstractSource> source, FetchSourceF fetch, PushSourceF push) {
    _source = std::move(source);
    _fetch = fetch;
    _push = push;
  }

  /// Outermost group, without the detached-subgroup check that root() applies
  ParametersD& outermost() {
    return _parent ? _parent->outermost() : *this;
  }

  /// Apply the highest-precedence input override (CLI/XML/...) for this field, if any
  void tryResolveOverride(TypeErasedFieldD& typeErasedField) {
    const OverrideSource::Query query{
      splitPrefix(_prefix), typeErasedField.name,
      typeErasedField.isNative, typeErasedField.isFlag};
    if (auto value = singleton::input().resolve(query)) {
      typeErasedField.parse(*this, *value);
    }
  }

  // Helper to get the typed field pointer, executing calculate if needed
  template <concepts::Field FIELD>
  FieldD<T,DESCRIPTOR,FIELD>* getFieldPtr() {
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    if (typeErasedField.calculate) {
      // Dependent parameters recompute on every access so they track their base
      if (typeErasedField.isBeingComputedRightNow) {
        throw std::runtime_error("Cyclic dependency detected in parameter: " + typeErasedField.name);
      }
      typeErasedField.isBeingComputedRightNow = true;
      try {
        typeErasedField.calculate(*this);
      } catch (...) {
        // Leaving the flag set would poison the field for every later access
        typeErasedField.isBeingComputedRightNow = false;
        throw;
      }
      typeErasedField.isBeingComputedRightNow = false;
      tryResolveOverride(typeErasedField);
    }
    else {
      if (!typeErasedField.field && !tryMaterializeFromSource<FIELD>()) {
        set<FIELD>(FIELD::template getInitialValue<T,DESCRIPTOR>());
      }
      if (!typeErasedField.resolved) {
        tryResolveOverride(typeErasedField);
        typeErasedField.resolved = true;
      }
    }
    return &static_cast<ConcreteFieldD<FIELD>*>(typeErasedField.field.get())->data;
  }

  explicit ParametersD(std::string prefix)
    : _prefix{std::move(prefix)}
  { }

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  ParametersD() = default;

  ParametersD(ParametersD&& rhs) noexcept
    : _map{std::move(rhs._map)}
    , _prefix{std::move(rhs._prefix)}
    , _children{std::move(rhs._children)}
    , _version{rhs._version}
    , _source{std::move(rhs._source)}
    , _fetch{rhs._fetch}
    , _push{rhs._push}
  {
    adoptChildren();
  }

  ParametersD& operator=(ParametersD&& rhs) noexcept {
    _map = std::move(rhs._map);
    _prefix = std::move(rhs._prefix);
    _children = std::move(rhs._children);
    _version = rhs._version;
    _source = std::move(rhs._source);
    _fetch = rhs._fetch;
    _push = rhs._push;
    adoptChildren();
    return *this;
  }

  ParametersD(const ParametersD& rhs)
    : _map{rhs._map}
    , _prefix{rhs._prefix}
    , _version{rhs._version}
    , _fetch{rhs._fetch}
    , _push{rhs._push}
  {
    for (const auto& [name, child] : rhs._children) {
      _children.emplace(name, std::make_unique<ParametersD>(*child));
    }
    adoptChildren();
    if (rhs._source) {
      attachSource(rhs._source->clone(), rhs._fetch, rhs._push);
    }
  }

  ParametersD& operator=(const ParametersD& rhs) {
    ParametersD copy{rhs};
    *this = std::move(copy);
    return *this;
  }

  /// Nested parameter subgroup, created on first use (--<prefix>.NAME, <Parameters><prefix>)
  ParametersD& operator[](const std::string& subgroup) {
    auto iter = _children.find(subgroup);
    if (iter == _children.end()) {
      std::string childPrefix = _prefix.empty() ? subgroup
                                                : _prefix + "." + subgroup;
      auto child = std::unique_ptr<ParametersD>(new ParametersD(std::move(childPrefix)));
      child->_parent = this;
      iter = _children.emplace(subgroup, std::move(child)).first;
    }
    return *iter->second;
  }

  /// Enclosing parameter group
  ParametersD& parent() {
    if (!_parent) {
      throw std::runtime_error("Parameter group '" + _prefix + "' has no parent");
    }
    return *_parent;
  }

  /// Outermost parameter group
  ParametersD& root() {
    if (_parent) {
      return _parent->root();
    }
    if (!_prefix.empty()) {
      throw std::runtime_error("Parameter group '" + _prefix + "' is detached from its tree");
    }
    return *this;
  }

  /// Bumped by every set<>; lets consumers cache extracted values
  std::size_t version() const {
    return _version;
  }

  /// Concrete values of FIELD_LIST, e.g. to hand an operator-style parameter set to a reporter
  template <typename FIELD_LIST>
  typename StaticParametersD<T,DESCRIPTOR>::template include<FIELD_LIST> extract() {
    typename StaticParametersD<T,DESCRIPTOR>::template include<FIELD_LIST> extracted;
    FIELD_LIST::for_each([&](auto id) {
      using FIELD = typename decltype(id)::type;
      extracted.template set<FIELD>(*this->getFieldPtr<FIELD>());
    });
    return extracted;
  }

  // Declare that we will use a parameter
  template <concepts::Field FIELD>
  void declare() {
    this->template get<FIELD>();
  }

  // Set parameter with a concrete value
  template <concepts::Field FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD> value) {
    if (!FIELD::template isValid<T,DESCRIPTOR,FIELD>(value)) {
      std::stringstream msg;
      msg << "Value " << value << " for " << fields::name<FIELD>() << " is invalid" << std::endl;
      throw std::invalid_argument(msg.str());
    }

    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];

    // If this set call is not triggered by recomputation it is a manual override
    if (!typeErasedField.isBeingComputedRightNow) {
      typeErasedField.calculate = nullptr;
    }

    typeErasedField.field = std::make_unique<ConcreteFieldD<FIELD>>(value);
    typeErasedField.fromSource = false;
    setupTypeErasedField<FIELD>(typeErasedField);
    ++_version;

    // The snapshot derives from its own values, so it needs this one too
    ParametersD& owner = outermost();
    if (owner._source && owner._push && !_materializing) {
      pushToSource<FIELD>(owner, value);
      owner.dropSourceValues();
    }
  }

  template <concepts::Field FIELD>
  requires (   DESCRIPTOR::template size<FIELD>() == 1
            && std::is_same_v<typename FIELD::template value_type<T>,std::string>)
  void set(const char* text) {
    set<FIELD>(std::string(text));
  }

  // Allows the conversion of double to ADf-type
  template <concepts::Field FIELD, typename U>
  requires (   DESCRIPTOR::template size<FIELD>() == 1
            && util::is_adf_v<typename FIELD::template value_type<T>>
            && std::same_as<std::remove_cvref_t<U>, BaseType<T>>)
  void set(U&& value) {
    set<FIELD>(FieldD<T,DESCRIPTOR,FIELD>(value));
  }

  /// Set dependent parameter; the lambda reads the set that evaluates it
  template <concepts::Field FIELD, typename F>
  requires (   std::invocable<F&,ParametersD&>
            && std::convertible_to<std::invoke_result_t<F&,ParametersD&>,
                                   FieldD<T,DESCRIPTOR,FIELD>>)
  void set(F&&) {
    static_assert(std::is_empty_v<std::decay_t<F>>,
      "a dependent parameter must not capture anything: read what it needs from the "
      "parameter set it is handed, so that copies and conversions evaluate correctly");
    TypeErasedFieldD& typeErasedField = _map[typeid(FIELD)];
    typeErasedField.field.reset();
    typeErasedField.calculate = [](ParametersD& self) {
      self.template set<FIELD>(FieldD<T,DESCRIPTOR,FIELD>(std::decay_t<F>{}(self)));
    };
    setupTypeErasedField<FIELD>(typeErasedField);
    ++_version;
  }

  template <concepts::Field FIELD>
  auto get() const {
    auto* ptr = const_cast<ParametersD*>(this)->getFieldPtr<FIELD>();
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      return (*ptr)[0];
    } else {
      return *ptr;
    }
  }


  template <concepts::Field FIELD>
  auto getOrFallback(const typename FieldD<T,DESCRIPTOR,FIELD>::value_type& fallback) const {
    if (provides<FIELD>()) {
      return get<FIELD>();
    } else {
      return fallback;
    }
  }

  template <concepts::Field FIELD>
  bool provides() const {
    auto it = _map.find(typeid(FIELD));
    if (it == _map.end()) {
      return false;
    }
    return it->second.field != nullptr || it->second.calculate != nullptr;
  }

  void print(std::ostream& clout = std::cout, const std::string& indent = "  ") {
    // Evaluating may insert entries, invalidating iterators but not pointers
    std::vector<std::pair<std::string,TypeErasedFieldD*>> entries;
    entries.reserve(_map.size());
    for (auto& [_, typeErasedField] : _map) {
      entries.emplace_back(typeErasedField.name, &typeErasedField);
    }
    std::sort(entries.begin(), entries.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    });

    io::text::KeyValueList list;
    for (auto& [name, typeErasedField] : entries) {
      list.item(name, typeErasedField->value(*this));
    }
    list.render(clout, indent);
    for (auto& [name, child] : _children) {
      clout << indent << name << ":" << std::endl;
      child->print(clout, indent + "  ");
    }
  }

  /// Same parameters in a different value type, each converted on first access
  /// Owns a snapshot, so dependent parameters arrive as values, not as derivations
  template <concepts::BaseType U>
  ParametersD<U,DESCRIPTOR> copyAs() const {
    using Target = ParametersD<U,DESCRIPTOR>;
    Target target;
    auto snapshot = std::make_unique<typename Target::template ConcreteSource<T>>(*this);
    target.attachSource(std::move(snapshot),
      [](typename Target::AbstractSource& source,
         const std::vector<std::string>& path,
         std::type_index field,
         void* out,
         const std::type_info& requested,
         std::size_t n) -> bool {
        auto& concrete = static_cast<typename Target::template ConcreteSource<T>&>(source);
        return concrete.data.exportAt(path, field, out, requested, n);
      },
      [](typename Target::AbstractSource& source,
         const std::vector<std::string>& path,
         std::type_index field,
         const void* in,
         const std::type_info& provided,
         std::size_t n) {
        auto& concrete = static_cast<typename Target::template ConcreteSource<T>&>(source);
        concrete.data.importAt(path, field, in, provided, n);
      });
    return target;
  }
};

namespace concepts {

/// Any (possibly cv-qualified or reference) ParametersD specialization
template <typename TYPE>
concept AnyParametersD = std::same_as<
  std::decay_t<TYPE>,
  olb::ParametersD<typename std::decay_t<TYPE>::value_t,
                   typename std::decay_t<TYPE>::descriptor_t>>;

}

}

#endif
