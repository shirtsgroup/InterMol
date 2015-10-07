echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

#if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
#    echo "This is a pull request. No deployment will be done."; exit 0
#fi
#
#
#if [[ "$TRAVIS_BRANCH" != "master" ]]; then
#    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
#fi


if [[ "2.7 3.4" =~ "$python" ]]; then
    binstar -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol $HOME/miniconda/conda-bld/linux-64/intermol-*
    conda convert $HOME/miniconda/conda-bld/linux-64/intermol-* -p all
    ls
    binstar -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol linux-32/intermol-*
    binstar -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol win-32/intermol-*
    binstar -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol win-64/intermol-*
    binstar -t "$BINSTAR_TOKEN"  upload --force --user shirtsgroup --package intermol osx-64/intermol-*
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi


