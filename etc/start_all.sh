#!/bin/sh -

echo downloading models:
/pombase-gocam-server/etc/get_gocams_by_taxon $mod_url /data

echo starting server:
exec /root/.local/bin/circusd --log-level debug /pombase-gocam-server/etc/circus.ini
