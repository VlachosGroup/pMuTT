# Conda build procedure

With Anaconda (or some other conda distribution) present in your environment, from this directory (`conda-pieces`) do:

```
$ mkdir /tmp/<user>_conda_bld
$ conda-build --croot=/tmp/<user>-conda-bld \
  --cache-dir=/tmp/<user>-conda-bld/tmp \
  --channel conda-forge \
  .
```

where `<user>` is your username on the machine in question.  The *conda-forge* channel is necessary to satisfy the dependency on ASE (ASE is not present on the default channels).

If all goes well, you should end up with a conda package at `/tmp/<user>_conda_bld/<os>/vlachos-group-mutt-1.0.0-py27_0.tar.bz2` (e.g. `<os>` might be `linux-64`).

## Regarding versioning

Note that the `meta.yaml` currently uses the `master` branch of the Git repository.  With proper version tagging happening in the repository, the `version` and `git_rev` parameters in `meta.yaml` should be in-sync (e.g. `version: 1.0.0` corresponding to `git_rev: v1.0.0`).


# Local conda channel

A local conda channel can be created to hold the resulting packages:

```
$ mkdir -p my_channel/{<os>,noarch}
$ cp /tmp/<user>_conda_bld/<os>/vlachos-group-mutt-1.0.0-py27_0.tar.bz2 my_channel/<os>
$ conda index my_channel/{<os>,noarch,}
```

Each time the package is replaced or a new one added, the directory should be re-indexed.

The channel is accessed by adding a `--channel` option to your `conda` commands:

```
$ conda search --channel=file:///<directory-containing-my_channel>/my_channel/ --override-channels
Loading channels: done
# Name                  Version           Build  Channel             
vlachos-group-mutt           1.0.0          py27_0  my_channel 
```

