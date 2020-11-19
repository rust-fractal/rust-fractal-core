mkdir -p output

for loc in 1e14 1e50 5e113 5e227 4e533 1e1086
do
    echo "${loc}.toml"
    echo "dim,real,user,sys,mem"
    for dim in 512 1024 2048 4096
    do
        /usr/bin/time --format="${dim},%e,%U,%S,%M" ../target/release/main -o "${dim}".toml "${loc}".toml 
    done
done