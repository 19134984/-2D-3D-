# This file is part of the OpenLB library
#
# Copyright (C) 2026 Shota Ito, Adrian Kummerlaender
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

import sys
import re

def format_string(s, INTROSPECTION=True):
    s = s.replace(" ", "")
    if INTROSPECTION:
        s = s.replace("double", "Expr")
        s = s.replace("float", "Expr")
    s = s.replace("olb::", "")
    return "".join(s)

if __name__ == "__main__":
    if len(sys.argv) > 2: # if pass a second argument disable introspection (double/float is not replaced by Expr)
        input_str = sys.argv[1]
        print(format_string(input_str, INTROSPECTION=False))
    elif len(sys.argv) > 1:
        input_str = sys.argv[1]
        print(format_string(input_str))
