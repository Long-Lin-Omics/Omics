

total=0
for f in $1/*.pod5; do
    n=$(pod5 inspect summary "$f" | grep Found |cut -d ' ' -f 4)
    total=$((total + n))
done
echo $total

