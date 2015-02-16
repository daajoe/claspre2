version="2.0.0"

to=release/claspre-$version

mkdir $to

cp -r app $to
cp -r libclasp $to
cp -r libprogram_opts $to
cp -r tools $to
cp -r CHANGES $to
cp -r configure.sh $to
cp -r COPYING $to
cp -r README $to

cd release
tar cvfz claspre-$version-src.tar.gz claspre-$version
cd ..
