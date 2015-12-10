#!/bin/bash
version=$1
#hadd -f diboson_$version.root WZ_inclusive_$version.root ZZ_inclusive_$version.root WpWp_$version.root WmWm_$version.root

hadd -f VVV_$version.root WWW_$version.root WWZ_$version.root WZZ_$version.root ZZZ_$version.root

hadd -f ttV_$version.root ttW_$version.root ttZ_$version.root ttWW_$version.root
