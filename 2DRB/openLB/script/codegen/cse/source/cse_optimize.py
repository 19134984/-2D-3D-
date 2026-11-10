# This file is part of the OpenLB library
#
# Copyright (C) 2024 Shota Ito, Adrian Kummerlaender
# E-mail contact: info@openlb.net
# The most recent release of OpenLB can be downloaded at
# <http://www.openlb.net/>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA  02110-1301, USA.

from sympy import *
from sympy.codegen.ast import CodeBlock, Assignment
from sympy.core.function import UndefinedFunction
from sympy.utilities.iterables import numbered_symbols
from cse_utils import SymbolGenerator, CodeBlockPrinter, code, cse
import re

_TEMP_NAME = re.compile(r'x\d+\Z')

def is_removable(assignment):
    # Only generated arithmetic temporaries (x0, x1, ...) may be dropped. cell[i]
    # and field-pointer writes are Symbols too but keep their structural names, and
    # computeRho/computeU/... calls may read or mutate the cell.
    has_call = any(isinstance(a.func, UndefinedFunction)
                   for a in assignment.rhs.atoms(Function))
    return bool(_TEMP_NAME.match(str(assignment.lhs))) and not has_call

def eliminate_dead_assignments(assignments, seed_live):
    # cse() never prunes its input, so shared subexpressions the final algebra does
    # not reference survive as dead temps. Mark outputs live, then keep only the
    # temporaries reachable from them. Names are preserved; dead lines just vanish.
    live = set(seed_live)
    for assignment in assignments:
        if not is_removable(assignment):
            live |= assignment.rhs.free_symbols
    keep = set()
    changed = True
    while changed:
        changed = False
        for assignment in assignments:
            if is_removable(assignment) and assignment.lhs in live and assignment.lhs not in keep:
                keep.add(assignment.lhs)
                live |= assignment.rhs.free_symbols
                changed = True
    kept = [a for a in assignments if not is_removable(a) or a.lhs in keep]
    return kept, live

def optimize_dynamics(inputFile, deterministic=False, target="Generic"):
    # Load expression tree from ouput file
    with open(inputFile, "r") as file:
        load_expressions = file.read()
    results = { }
    exec(load_expressions, globals(), results)

    # Operator information
    collisionO = results['collisionO']
    # Assignments
    cell_assignments = results['cell_assignments']
    return_assignments = results['return_assignments']
    fields_assignments = results['fields_assignments']
    # Symbols
    params_symbols = results['params_symbols']
    shared_subexpressions = results.get('cse_subexpressions', [ ])

    # Retrieve dynamics tuple information
    collisionO = collisionO.replace("olb::", "")
    collisionO = collisionO.replace("Expr", "T")
    collisionO = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', collisionO)

    # Current limitations of sympy assignments
    generator = iter(SymbolGenerator('x'))

    # Rename the shared subexpressions onto the x-generator (so they are declared
    # auto) and rewrite later references; bindings arrive in dependency order.
    shared_assignments = [ ]
    shared_rename = { }
    for assignment in shared_subexpressions:
        alias = next(generator)
        shared_assignments.append(Assignment(alias, assignment.rhs.xreplace(shared_rename)))
        shared_rename[assignment.lhs] = alias

    if shared_rename:
        cell_assignments   = [Assignment(a.lhs, a.rhs.xreplace(shared_rename)) for a in cell_assignments]
        fields_assignments = [Assignment(a.lhs, a.rhs.xreplace(shared_rename)) for a in fields_assignments]
        return_assignments = [Assignment(a.lhs, a.rhs.xreplace(shared_rename)) for a in return_assignments]

    # Combine all assignments
    assignments = cell_assignments
    if deterministic:
        assignments.sort(key=lambda x: str(x.lhs))
    optional_symbols = params_symbols
    if deterministic:
        optional_symbols = sorted(optional_symbols, key=lambda x: str(x))

    cell_backsubstitutions = [ ]
    sub_assignments = [ ]
    for assignment in assignments:
        alias = next(generator)
        cell_backsubstitutions.append((alias, assignment.lhs))
        sub_assignments.append(Assignment(alias, assignment.rhs))

    # Creating aliases and assignments for parameter calls
    optional_assignments = [ ]
    optional_substitutions = [ ]
    for symbol in optional_symbols:
        alias = next(generator)
        optional_assignments.append(Assignment(alias, symbol))
        optional_substitutions.append((symbol, alias))

    # Create code block to perform cse
    block = CodeBlock(*shared_assignments, *sub_assignments, *fields_assignments, *return_assignments)
    block = block.subs(optional_substitutions)

    # Perform cse
    block = cse(block, generator, deterministic=deterministic)

    # These are required due to current CSE limititations
    substitutions = [ ]
    post_collision = [ ]
    simple_post_collision = True
    for assignment in block.args:
        if assignment.lhs in cell_backsubstitutions[0]:
            atoms = assignment.rhs.free_symbols
            for (lhs, rhs) in cell_backsubstitutions:
                if lhs in atoms or rhs in atoms:
                    simple_post_collision = False
    if simple_post_collision:
        substitutions = cell_backsubstitutions
    else:
        for subs in cell_backsubstitutions:
            post_collision.append(Assignment(subs[1], subs[0]))

    # Apply back-substitution of placeholders
    block = block.subs(substitutions)

    # Detect which optional assignments are required to close the symbolic collision
    optional = [ ]
    free_syms = block.free_symbols
    if deterministic:
        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower()
                    for text in re.split('([0-9]+)', str(s))]
        free_syms = sorted(free_syms, key=natural_sort_key)

    for symbol in free_syms:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    optional.append(assgn)

    # Drop shared subexpressions the final collision never references (dead temps).
    body = list(block.args[:-2])
    seed_live = set(block.args[-2].rhs.free_symbols) | set(block.args[-1].rhs.free_symbols)
    for assignment in post_collision:
        seed_live |= assignment.rhs.free_symbols
    body, live = eliminate_dead_assignments(body, seed_live)
    optional = [assignment for assignment in optional if assignment.lhs in live]

    # return cse optimized expression
    return '\n'.join([
        "template <typename T, typename... FIELDS>",
        f"""struct CSE<{ collisionO }, olb::CodegenTarget::{ target }> {{""",
        "template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>",
        "CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {",
        code(CodeBlock(*optional, *body, *post_collision), symbols=generator.symbols),
        "return { %s, %s };" % (code(block.args[-2].rhs), code(block.args[-1].rhs)),
        "}\n};"
    ])

def optimize_operator(inputFile):
    # Load expression tree from ouput file
    with open(inputFile, "r") as file:
        load_expressions = file.read()
    results = { }
    exec(load_expressions, globals(), results)

    # Operator information
    operatorO = results['operatorO']
    descriptor = results['descriptor']
    isOperatorWithParameter = results['isOperatorWithParameter']
    # Assignments
    cell_assignments = results['cell_assignments']
    fields_assignments = results['fields_assignments']
    function_calls = results['function_calls']
    # Symbols
    params_symbols = results['params_symbols']
    shared_subexpressions = results.get('cse_subexpressions', [ ])

    operatorO = operatorO.replace("olb::", "")
    operatorDependsOnT = "Expr" in operatorO
    operatorO = operatorO.replace("Expr", "T")
    operatorO = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', operatorO)
    descriptor = re.sub(r'(descriptors::D\dQ\d{1,2})<.*?>', r'\1<FIELDS...>', descriptor)

    # Current limitations of sympy assignments
    generator = iter(SymbolGenerator('x'))

    # Rename the shared subexpressions onto the x-generator; they reference symbols
    # from function_calls, so they are placed after them in the code block below.
    shared_assignments = [ ]
    shared_rename = { }
    for assignment in shared_subexpressions:
        alias = next(generator)
        shared_assignments.append(Assignment(alias, assignment.rhs.xreplace(shared_rename)))
        shared_rename[assignment.lhs] = alias

    if shared_rename:
        cell_assignments   = [Assignment(a.lhs, a.rhs.xreplace(shared_rename)) for a in cell_assignments]
        fields_assignments = [Assignment(a.lhs, a.rhs.xreplace(shared_rename)) for a in fields_assignments]

    # Combine all assignments
    assignments = cell_assignments
    optional_symbols = params_symbols

    cell_backsubstitutions = [ ]
    sub_assignments = [ ]
    for assignment in assignments:
        alias = next(generator)
        cell_backsubstitutions.append((alias, assignment.lhs))
        sub_assignments.append(Assignment(alias, assignment.rhs))

    # Creating aliases and assignments for parameter calls
    optional_assignments = [ ]
    optional_substitutions = [ ]
    for symbol in optional_symbols:
        alias = next(generator)
        optional_assignments.append(Assignment(alias, symbol))
        optional_substitutions.append((symbol, alias))

    # Create code block to perform cse
    block = CodeBlock(*function_calls, *shared_assignments, *sub_assignments, *fields_assignments)
    block = block.subs(optional_substitutions)

    # Perform cse
    block = cse(block, generator)

    # These are required due to current CSE limititations
    substitutions = [ ]
    post_apply = [ ]
    simple_post_apply = True
    for assignment in block.args:
        if assignment.lhs in cell_backsubstitutions[0]:
            atoms = assignment.rhs.free_symbols
            for (lhs, rhs) in cell_backsubstitutions:
                if lhs in atoms or rhs in atoms:
                    simple_post_apply = False
    if simple_post_apply:
        substitutions = cell_backsubstitutions
    else:
        for subs in cell_backsubstitutions:
            post_apply.append(Assignment(subs[1], subs[0]))

    # Apply back-substitution of placeholders
    block = block.subs(substitutions)

    # Detect which optional assignments are required to close the symbolic apply
    optional = [ ]
    for symbol in block.free_symbols:
        if symbol in { assgn.lhs for assgn in optional_assignments }:
            for assgn in optional_assignments:
                if assgn.lhs is symbol:
                    optional.append(assgn)

    # Remove helper assignments to ensure symbol availability
    args = [ ]
    for assignment in block.args:
        if hasattr(assignment.rhs, 'name'):
            if "EXPOSE_ARRAY_ELEMENT" not in assignment.rhs.name:
                args.append(assignment)
        else:
            args.append(assignment)

    # Drop shared subexpressions the final apply never references (dead temps).
    seed_live = set()
    for assignment in post_apply:
        seed_live |= assignment.rhs.free_symbols
    kept, live = eliminate_dead_assignments(args, seed_live)
    block = CodeBlock(*kept)
    optional = [assignment for assignment in optional if assignment.lhs in live]

    # return cse optimized expression
    header = [ ]
    if isOperatorWithParameter:
        header = "\n".join([
            "template <concepts::Cell CELL, concepts::Parameters PARAMETERS>",
            "void apply(CELL& cell, PARAMETERS& parameters) any_platform {"
        ])
    else:
        header = "\n".join([
            "template <concepts::Cell CELL>",
            "void apply(CELL& cell) any_platform {"
        ])

    template_header = "template <typename... FIELDS>"
    if operatorDependsOnT:
        template_header = "template <typename T, typename... FIELDS>"
    return '\n'.join([
        template_header,
        f"""struct CSE_O<{ operatorO },{ descriptor }> {{""",
        header,
        "using V = typename CELL::value_t;",
        "using DESCRIPTOR = typename CELL::descriptor_t;",
        code(CodeBlock(*optional, *block.args, *post_apply), symbols=generator.symbols),
        "}\n};"
    ])
