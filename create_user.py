#!/usr/bin/env python3

# Copyright (C) 2019 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import argparse


def create_LDIF_entry(cli_args):
    print(cli_args.name)
    print(cli_args.familiy_name)

parser = argparse.ArgumentParser(prog="create_user", formatter_class=argparse.RawDescriptionHelpFormatter,
                                 fromfile_prefix_chars="@",
                                 description="Contact: tobias.jakobi@med.uni-heidelberg.de")

group = parser.add_argument_group("Input")

group.add_argument("-p",
                   "--password",
                   dest="password",
                   help="Initial user password"
                   )

group.add_argument("-U",
                   "--user-id",
                   dest="user_id",
                   help="Numerical user ID to start with (for batch creations)"
                   )

group.add_argument("-g",
                   "--group",
                   dest="group_id",
                   help="Numerical group to add the user to",
                   required=True
                   )

group.add_argument("-e",
                   "--email",
                   dest="email",
                   help="User email",
                   required=True
                   )

group.add_argument("-n",
                   "--name",
                   dest="name",
                   help="Name of the user",
                   required=True
                   )

group.add_argument("-f",
                   "--familiy-name",
                   dest="family_name",
                   help="Familiy name of the user",
                   required=True
                   )

args = parser.parse_args()

create_LDIF_entry(args)
