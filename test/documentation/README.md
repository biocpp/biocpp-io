# bio documentation

Currently, we can only build the documentation on *nix systems.

We offer two versions of our documentation, one intended for the user (doc_usr) and one intended for the
library-developer and maintainer (doc_dev) of bio, which contains the documentation of internals.

## How to configure:

```bash
mkdir <bio-build-dir>/documentation
cd <bio-build-dir>/documentation

cmake <bio-dir>/test/documentation
```

## How to build:

Prerequisites: configuring the documentation

```bash
cd <bio-build-dir>/documentation

# build user and developer documentation
cmake --build .

# or build user documentation only
cmake --build . --target doc_usr

# or build developer documentation only
cmake --build . --target doc_dev
```

## How to test:

Prerequisites: building the documentation

```bash
cd <bio-build-dir>/documentation

# test user and developer documentation
ctest --output-on-failure --progress

# or test user documentation only
ctest -R doc_usr --output-on-failure --progress

# or test developer documentation only
ctest -R doc_dev --output-on-failure --progress
```
