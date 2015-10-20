echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "develop" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi


if [[ "2.7 3.4" =~ "$python" ]]; then
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol-dev $HOME/miniconda/conda-bld/linux-64/intermol-*
    conda convert $HOME/miniconda/conda-bld/linux-64/intermol-* -p all
    pwd
    ls *-64/*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol-dev linux-32/intermol-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol-dev win-32/intermol-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol-dev win-64/intermol-*
    anaconda -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol-dev osx-64/intermol-*
fi
