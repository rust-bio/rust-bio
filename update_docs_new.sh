rm -r doc
cargo doc --no-deps
cp -r target/doc doc
git checkout gh-pages
git add -A doc
git commit -a
git push
git checkout master
