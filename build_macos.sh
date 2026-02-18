#!/bin/bash

# Exit on error
set -e

echo "🚀 Starting build process for fasterx..."

# 0. Clean old build artifacts to prevent recursive bundling
echo "🧹 Cleaning old build artifacts..."
rm -rf dist build

# 1. Bundle the Dioxus app in release mode
echo "📦 Bundling Dioxus app..."
export MACOSX_DEPLOYMENT_TARGET=12.0
dx bundle --release --platform desktop

# The output path for Dioxus 0.5/0.6 desktop apps is usually in dist/
APP_PATH="dist/bundle/macos/fasterx.app"

# If the path above doesn't exist, we'll try to find it
if [ ! -d "$APP_PATH" ]; then
    echo "🔍 Searching for the .app bundle..."
    FOUND_PATH=$(find dist -name "*.app" -type d | head -n 1)
    if [ -n "$FOUND_PATH" ]; then
        APP_PATH="$FOUND_PATH"
        echo "✅ Found app at: $APP_PATH"
    else
        echo "❌ Could not find the .app bundle. Please check 'dx build' output."
        exit 1
    fi
fi

# 2. Stripping security attributes
echo "🔐 Stripping security attributes..."

# Remove quarantine attribute (Gatekeeper bypass)
echo "   - Removing quarantine attributes..."
xattr -cr "$APP_PATH"

# Ad-hoc sign the bundle
# This is necessary for ARM64 (Apple Silicon) to run the binary at all
echo "   - Applying ad-hoc signature..."
codesign --force --deep --sign - "$APP_PATH"

# 3. Setting specific version metadata
echo "📝 Setting version metadata in Info.plist..."
plutil -replace LSMinimumSystemVersion -string "12.0" "$APP_PATH/Contents/Info.plist"

echo "✨ Build complete! The executable is ready at: $APP_PATH"
echo "👉 You can now share this .app bundle. Users might still need to right-click 'Open' the first time."
