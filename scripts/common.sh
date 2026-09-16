# Shared shell helpers, sourced by every run script.

# `mkdir -p /eos/...` fails on the EOS FUSE mount when an existing ancestor is only traversable via
# EOS ACLs (coreutils tries chdir/mkdir on each ancestor). Create only the components that are missing.
mkdirp() {
  local d="$1" todo=()
  while [ ! -d "$d" ]; do todo=("$d" "${todo[@]}"); d="$(dirname "$d")"; done
  local x; for x in "${todo[@]}"; do mkdir "$x"; done
}

# Bind mounts every containerised step needs on this site: CVMFS (software), EOS + the Kerberos
# credential cache (EOS FUSE authorises container processes through the host's krb5 ticket in
# /run/user/<uid>; without that bind, EOS is read/write-denied inside apptainer), and /tmp.
apptainer_binds() {
  local b="-B /cvmfs -B /tmp"
  [ -d /eos ] && b="$b -B /eos"
  [ -d "/run/user/$(id -u)" ] && b="$b -B /run/user/$(id -u)"
  echo "$b"
}
