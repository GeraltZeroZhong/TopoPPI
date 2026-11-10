#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: $0 VERSION ARCHITECTURE ENV_PREFIX OPTCUTS_BIN OUTPUT_DIR" >&2
    exit 2
fi

version="$1"
architecture="$2"
environment_prefix="$3"
optcuts_bin="$4"
output_dir="$5"

case "$architecture" in
    arm64|x86_64) ;;
    *) echo "Unsupported macOS architecture: $architecture" >&2; exit 2 ;;
esac

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
build_root="$(mktemp -d "${TMPDIR:-/tmp}/topoppi-macos-app.XXXXXX")"
app_bundle="$build_root/TopoPPI.app"
contents_dir="$app_bundle/Contents"
resources_dir="$contents_dir/Resources"
iconset_dir="$build_root/topoppi.iconset"
disk_image="TopoPPI-${version}-macos-${architecture}.dmg"

cleanup() {
    rm -rf -- "$build_root"
}
trap cleanup EXIT

mkdir -p "$contents_dir/MacOS" "$resources_dir" "$iconset_dir" "$output_dir"
install -m 755 "$script_dir/TopoPPI" "$contents_dir/MacOS/launch.sh"
xcrun clang \
    -Os \
    -mmacosx-version-min=12.0 \
    "$script_dir/launcher.c" \
    -o "$contents_dir/MacOS/TopoPPI"
test "$(lipo -archs "$contents_dir/MacOS/TopoPPI")" = "$architecture"
sed \
    -e "s/@VERSION@/$version/g" \
    -e "s/@ARCHITECTURE@/$architecture/g" \
    "$script_dir/Info.plist.in" > "$contents_dir/Info.plist"

icon_source="$repo_root/src/topoppi/assets/topoppi.png"
for icon_spec in \
    "16 icon_16x16.png" \
    "32 icon_16x16@2x.png" \
    "32 icon_32x32.png" \
    "64 icon_32x32@2x.png" \
    "128 icon_128x128.png" \
    "256 icon_128x128@2x.png" \
    "256 icon_256x256.png" \
    "512 icon_256x256@2x.png" \
    "512 icon_512x512.png" \
    "1024 icon_512x512@2x.png"
do
    read -r pixels filename <<< "$icon_spec"
    sips -z "$pixels" "$pixels" "$icon_source" --out "$iconset_dir/$filename" >/dev/null
done
iconutil -c icns "$iconset_dir" -o "$resources_dir/topoppi.icns"

"$environment_prefix/bin/python" "$script_dir/collect_licenses.py" "$environment_prefix" "$resources_dir"
install -m 644 "$repo_root/LICENSE" "$resources_dir/TopoPPI-LICENSE.txt"
install -m 644 "$repo_root/tools/OptCuts/LICENSE.txt" "$resources_dir/OptCuts-LICENSE.txt"
install -m 644 "$repo_root/tools/OptCuts/NOTICE.md" "$resources_dir/OptCuts-NOTICE.md"
install -m 644 "$repo_root/tools/OptCuts/THIRD_PARTY_LICENSES.txt" "$resources_dir/OptCuts-THIRD-PARTY-LICENSES.txt"

raw_environment="$build_root/environment.raw.tar.gz"
packed_environment="$build_root/environment"
conda-pack \
    -p "$environment_prefix" \
    -o "$raw_environment" \
    --exclude 'lib/python*/site-packages/topoppi-*.dist-info/direct_url.json' \
    --exclude 'lib/python*/site-packages/topoppi_affinity-*.dist-info/direct_url.json' \
    --exclude 'lib/python*/site-packages/topoppi_hotspot-*.dist-info/direct_url.json'
mkdir -p "$packed_environment"
tar -xzf "$raw_environment" -C "$packed_environment"
find "$packed_environment" -type f \( -name '*.pyc' -o -name '*.pyo' \) -delete
find "$packed_environment" -type d -name __pycache__ -prune -exec rm -rf {} +
rm -rf "$packed_environment/conda-meta"

own_runtime_paths=()
while IFS= read -r -d '' path; do
    own_runtime_paths+=("$path")
done < <(
    find "$packed_environment" -type d \
        \( -path '*/site-packages/topoppi' -o -path '*/site-packages/topoppi_chemistry' \
           -o -path '*/site-packages/topoppi_representation' \
           -o -path '*/site-packages/topoppi_affinity' -o -path '*/site-packages/topoppi_affinity_tools' \
           -o -path '*/site-packages/topoppi_hotspot' \
           -o -path '*/site-packages/topoppi-*.dist-info' \
           -o -path '*/site-packages/topoppi_affinity-*.dist-info' \
           -o -path '*/site-packages/topoppi_hotspot-*.dist-info' \) \
        -print0
)
private_paths="$({
    LC_ALL=C grep -RIl -a -F \
        -e "$repo_root" \
        "${own_runtime_paths[@]}" \
        "$contents_dir/MacOS" \
        "$optcuts_bin" || true
} | sort -u)"
if [[ -n "$private_paths" ]]; then
    echo "Private build paths remain in the bundled runtime:" >&2
    echo "$private_paths" >&2
    exit 1
fi
tar -czf "$resources_dir/environment.tar.gz" -C "$packed_environment" .
install -m 755 "$optcuts_bin" "$resources_dir/OptCuts_bin"

plutil -lint "$contents_dir/Info.plist"
codesign --force --deep --sign - "$app_bundle"
codesign --verify --deep --strict "$app_bundle"

disk_root="$build_root/disk"
mkdir -p "$disk_root"
ditto "$app_bundle" "$disk_root/TopoPPI.app"
ln -s /Applications "$disk_root/Applications"
hdiutil create \
    -volname "TopoPPI $version" \
    -srcfolder "$disk_root" \
    -format UDZO \
    -ov \
    "$output_dir/$disk_image"

echo "$output_dir/$disk_image"
