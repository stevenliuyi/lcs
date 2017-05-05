#!/bin/sh

# clone gh-pages branch
git clone -b gh-pages https://github.com/stevenliuyi/lcs code_docs

# clear gh-pages branch
rm -rf code_docs/*

# generate documentation
doxygen docs/Doxyfile

# move documentation files to gh-pages branch
mv docs/html/* code_docs/

cd code_docs

# configure git
git config user.name "stevenliuyi"
git config user.email "stevenliuyi@gmail.com"

# add documentations and push back to gh-pages
git add --all
git commit -m "Doxygen deployment - Travis build: ${TRAVIS_BUILD_NUMBER}"
git push https://${GITHUB_TOKEN}@github.com/stevenliuyi/lcs gh-pages:gh-pages --force
