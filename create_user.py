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
import codecs
import hashlib
import os
import random
import re
import textwrap


# Generate XKCD style passwords
# from: https://stackoverflow.com/a/9368832/1900920

# apt-get install wbritish
def random_words(num, dictionary="/usr/share/dict/british-english"):
    r = random.SystemRandom()  # i.e. preferably not pseudo-random
    f = open(dictionary, "r")
    count = 0
    chosen = []
    for i in range(num):
        chosen.append("")
    prog = re.compile("^[a-z]{5,9}$")  # reasonable length, no proper nouns
    if f:
        for word in f:
            if prog.match(word):
                for i in range(num):  # generate all words in one pass thru file
                    if r.randint(0, count) == 0:
                        chosen[i] = word.strip()
                count += 1
    return chosen


def gen_password(num=1):
    return "".join(random_words(num))


# end password generation


# Get SSHA value
# from:  https://www.openldap.org/faq/data/cache/347.html:

from base64 import urlsafe_b64encode as encode
from base64 import urlsafe_b64decode as decode


def get_ssha_password(password):
    salt = os.urandom(4)
    h = hashlib.sha1(password.encode('utf-8'))
    h.update(salt)
    return "{SSHA}" + str(encode(h.digest() + salt))


def check_password(challenge_password, password):
    challenge_bytes = decode(challenge_password[6:])
    digest = challenge_bytes[:20]
    salt = challenge_bytes[20:]
    hr = hashlib.sha1(password)
    hr.update(salt)
    return digest == hr.digest()


# End SSHA functions


def username_from_full_name(cli_args):

    if not cli_args.username:

        name = cli_args.name
        family_name = cli_args.family_name
        return (name[0] + family_name).lower()

    else:
        return cli_args.username


def samba_ntpassword_from_password(password):
    nt_password = hashlib.new('md4', password.encode('utf-16le')).digest()
    nt_password = codecs.encode(nt_password, 'hex_codec').decode('utf-8').upper()
    return nt_password


def samba_uid_from_uid(uid, domain):
    return domain + "-" + str(uid * 2 + 1000)


def samba_guid_from_group_id(gid, domain):
    return domain + "-" + str(gid * 2 + 1001)


def generate_ldif_file(cli_args):
    if not cli_args.password:
        pw = gen_password()
        cli_args.password = pw
        print("# No password given, auto-generated: " + str(pw))
    else:
        print("# Given password: " + cli_args.password)

    print("version: 1")

    print("dn: cn=" + cli_args.name + " " + cli_args.family_name + ",ou=Users,dc=dieterichlab,dc=org")
    print("cn: " + cli_args.name + " " + cli_args.family_name)
    print("displayname: " + cli_args.name + " " + cli_args.family_name)
    print("gidnumber: " + str(cli_args.group_id))
    print("givenname: " + cli_args.name)
    print("homedirectory: /home/" + username_from_full_name(cli_args))
    print("loginshell: /bin/bash")
    print("mail: " + cli_args.email)

    print(textwrap.dedent("""\
    objectclass: shadowAccount
    objectclass: sambaSamAccount
    objectclass: posixAccount
    objectclass: inetOrgPerson
    objectclass: organizationalPerson
    objectclass: person
    """))

    print("sambaacctflags: [XU         ]")
    print("sambadomainname: " + cli_args.domain)
    print("sambahomedrive: U:")
    print("sambantpassword: " + samba_ntpassword_from_password(cli_args.password))
    print("sambaprimarygroupsid: " + samba_guid_from_group_id(cli_args.group_id, cli_args.samba_id))
    print("sambapwdlastset: 1554477616")
    print("sambasid: " + samba_uid_from_uid(cli_args.user_id, cli_args.samba_id))

    print(textwrap.dedent("""\
    shadowinactive: 10
    shadowlastchange: 17991
    shadowmax: 365
    shadowmin: 1
    shadowwarning: 10
    """))
    print("sn: " + cli_args.family_name)
    print("uid: " + username_from_full_name(cli_args))
    print("uidnumber: " + str(cli_args.user_id))
    print("userpassword: " + get_ssha_password(cli_args.password))


def create_ldif_entry(cli_args):

    if cli_args.domain == "AZ3":
        cli_args.samba_id = "S-1-5-21-1426298215-3934214462-2063419693"
    else:
        cli_args.samba_id = "S-1-5-21-3632671680-3219116354-420167436"

    cli_args.name = cli_args.name.rstrip()
    cli_args.family_name = cli_args.family_name.rstrip()

    generate_ldif_file(cli_args)

# main script starts here


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
                   help="Numerical user ID to start with (for batch creations)",
                   required=True,
                   type=int
                   )

group.add_argument("-u",
                   "--user-name",
                   dest="username",
                   help="UNIX user name (overides auto-generated username)"
                   )

group.add_argument("-g",
                   "--group",
                   dest="group_id",
                   help="Numerical group to add the user to",
                   required=True,
                   type=int
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
                   "--family-name",
                   dest="family_name",
                   help="Family name of the user",
                   required=True
                   )

group.add_argument("-d",
                   "--domain",
                   dest="domain",
                   help="Dieterichlab or AZ3 domain?",
                   choices=("Dieterichlab", "AZ3"),
                   required=True
                   )
args = parser.parse_args()

create_ldif_entry(args)
