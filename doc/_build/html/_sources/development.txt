============
Development
============

Contributing
------------
Contributions are welcome and they are greatly appreciated! Every little bit
helps and credit will always be given.

Git Flow
--------
Because we are supporting multiple molecular dynamics engines that should all
work independently, we try to keep development of each engine in a separate
branch until the basics are working.

To this end, we've started working with the `git flow branching model
<http://nvie.com/posts/a-successful-git-branching-model/>`_. The basic things
to know are:

1. The master branch is strictly used for releases.
2. The develop branch is where (!!!) development happens.
3. When we start working on a new engine, we create a feature branch. E.g.,
   at the time of this writing, there are branches called feature/lammps and
   feature/desmond. Once the overall structure in this branch is fairly
   stable and has a good amounts of tests, we merge it into develop.

So what do you, the interested developer, need to know?

1. Don't make pull requests against master.
2. Choose either develop or the appropriate feature branch to pull against.

For more reading and a neat tool to help with branching see:

http://jeffkreeftmeijer.com/2010/why-arent-you-using-git-flow/

https://github.com/nvie/gitflow

http://danielkummer.github.io/git-flow-cheatsheet/
