#!/bin/bash
# https://raw.githubusercontent.com/lotharschulz/sphinx-pages/master/sphinx-2-gh-pages-converter.sh

set -euo pipefail
IFS=$'\n\t'

buildDirectory=docs/_build

# get a clean master branch assuming
git checkout master
git pull origin master
git clean -df
git checkout -- .
git fetch --all

# build html docs from sphinx files
sphinx-build -b html . "$buildDirectory"

# create or use orphaned gh-pages branch
branch_name=gh-pages
if [ $(git branch --list "$branch_name") ]
then
	git stash
	git checkout $branch_name
	git pull origin $branch_name
	#git stash apply
	git checkout stash -- . # force git stash to overwrite added files
else
	git checkout --orphan "$branch_name"
fi

if [ -d "$buildDirectory" ]
then
	ls | grep -v _build | xargs rm -r
	mv _build/* . && rm -rf _build
	git add .
	git commit -m "new pages version $(date)"
	git push origin gh-pages
	# github.com recognizes gh-pages branch and create pages
	# url scheme https//:[github-handle].github.io/[repository]
else
	echo "directory $buildDirectory does not exists"
fi
