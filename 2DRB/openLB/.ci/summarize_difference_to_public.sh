#!/usr/bin/env bash

set -eo pipefail

curl -fsS "https://gitlab.com/openlb/release/-/archive/master/release-master.tar.gz" --output olb-release.tgz

tar xf olb-release.tgz
tar xf olb-master.tgz

pushd release-master
git init -q
git config user.name "OpenLB CI Worker"
git config user.email "ci@openlb.net"
git add -Af .
git commit -q -m "Old release"
git rm -rq --ignore-unmatch .
popd

cp -a olb-master/. release-master/

pushd release-master
git add -Af .
git diff --binary --staged --full-index > ../update-release-to-master.patch
git diff --staged --stat -- 'src/cse' > ../omitted-from-summary.txt
git diff --staged -- . ':(exclude)src/cse' > ../summary.patch
popd

if [ -s summary.patch ]; then
    diff2html --summary open --diffMaxChanges 400 --diffMaxLineLength 500 --title "Changes since last commit to public repository (generated src/cse omitted, binaries listed only)" --file summary.html --input file -- summary.patch
else
    echo "<html><head/><body><h1>Nothing changed outside src/cse</h1></body></html>" > summary.html
fi

echo "omitted from summary, see omitted-from-summary.txt:"
tail -1 omitted-from-summary.txt
