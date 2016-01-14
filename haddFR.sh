#!/bin/bash
version=$1
#hadd -f diboson_$version.root WZ_inclusive_$version.root ZZ_inclusive_$version.root WpWp_$version.root WmWm_$version.root

hadd -f VVV$version.root WWW$version.root WWZ$version.root WZZ$version.root ZZZ$version.root

hadd -f ttV$version.root ttW$version.root ttZ$version.root ttWW$version.root
