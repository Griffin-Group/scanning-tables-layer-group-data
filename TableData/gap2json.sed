#!/usr/bin/sed -E -f
s/tables := //
s/ := / : /g
s/ rec\(/ {/g
s/\)/\}/g
s/'/"/g
s/"?([a-zA-Z_]+)"?/"\1"/g
s/"?fail"?/null/g
s/"true"/true/g
s/"false"/false/g
s/;//g
s#(-?[1-9]+/[1-9]+)#"\1"#g
