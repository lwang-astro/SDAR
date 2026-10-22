version=`git rev-list --count HEAD`
version=`expr $version + 1`
# Preserve the 'e' experiment suffix across bumps. Appending it by hand
# after the script was error-prone (2026-10-22: a manual overwrite produced
# 469f instead of 470e; the script itself was never buggy).
suffix=''
[ -f VERSION ] && grep -q 'e$' VERSION && suffix='e'
echo "${version}${suffix}" >VERSION
