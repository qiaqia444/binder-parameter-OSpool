#!/bin/bash

# Extract Results Management Script
# Extracts tar.gz files from downloaded_results/ to appropriate directories

set -e

DOWNLOADED_DIR="downloaded_results"
OUTPUT_DIR="output"
LOGS_DIR="logs"

echo "============================================================"
echo "EXTRACT RESULTS MANAGEMENT"
echo "============================================================"

# Function to extract results
extract_archive() {
    local archive_name="$1"
    local archive_path="${DOWNLOADED_DIR}/${archive_name}"
    
    if [[ ! -f "$archive_path" ]]; then
        echo "❌ Archive not found: $archive_path"
        return 1
    fi
    
    echo "📦 Extracting: $archive_name"
    
    # Create temporary extraction directory
    local temp_dir="temp_extract_$$"
    mkdir -p "$temp_dir"
    
    # Extract to temporary directory
    tar -xzf "$archive_path" -C "$temp_dir"
    
    # Move extracted files to appropriate directories
    if [[ -d "$temp_dir/jobs/output" ]]; then
        echo "  → Moving simulation results to $OUTPUT_DIR/"
        mkdir -p "$OUTPUT_DIR"
        cp -r "$temp_dir/jobs/output/"* "$OUTPUT_DIR/" 2>/dev/null || true
    fi
    
    if [[ -d "$temp_dir/output" ]]; then
        echo "  → Moving simulation results to $OUTPUT_DIR/"
        mkdir -p "$OUTPUT_DIR"
        cp -r "$temp_dir/output/"* "$OUTPUT_DIR/" 2>/dev/null || true
    fi
    
    if [[ -d "$temp_dir/jobs/logs" ]]; then
        echo "  → Moving logs to $LOGS_DIR/"
        mkdir -p "$LOGS_DIR"
        cp -r "$temp_dir/jobs/logs/"* "$LOGS_DIR/" 2>/dev/null || true
    fi
    
    if [[ -d "$temp_dir/logs" ]]; then
        echo "  → Moving logs to $LOGS_DIR/"
        mkdir -p "$LOGS_DIR"
        cp -r "$temp_dir/logs/"* "$LOGS_DIR/" 2>/dev/null || true
    fi
    
    # Clean up temporary directory
    rm -rf "$temp_dir"
    
    echo "  ✅ Extraction complete"
}

# Function to list available archives
list_archives() {
    echo "Available archives in $DOWNLOADED_DIR/:"
    echo "----------------------------------------"
    for archive in "$DOWNLOADED_DIR"/*.tar.gz; do
        if [[ -f "$archive" ]]; then
            local basename=$(basename "$archive")
            local size=$(du -h "$archive" | cut -f1)
            echo "  📦 $basename ($size)"
        fi
    done
}

# Function to show current status
show_status() {
    echo "Current Status:"
    echo "---------------"
    echo "📁 Downloaded archives: $(ls -1 "$DOWNLOADED_DIR"/*.tar.gz 2>/dev/null | wc -l)"
    echo "📄 JSON results: $(ls -1 "$OUTPUT_DIR"/*.json 2>/dev/null | wc -l)"
    echo "📋 Log files: $(ls -1 "$LOGS_DIR"/* 2>/dev/null | wc -l)"
}

# Main menu
case "${1:-}" in
    "list")
        list_archives
        ;;
    "extract")
        if [[ -z "$2" ]]; then
            echo "Usage: $0 extract <archive_name>"
            echo "Available archives:"
            list_archives
            exit 1
        fi
        extract_archive "$2"
        ;;
    "extract-all")
        echo "Extracting all archives..."
        for archive in "$DOWNLOADED_DIR"/*.tar.gz; do
            if [[ -f "$archive" ]]; then
                extract_archive "$(basename "$archive")"
            fi
        done
        ;;
    "status")
        show_status
        ;;
    "help"|*)
        echo "Extract Results Management Script"
        echo ""
        echo "Usage:"
        echo "  $0 list           - List available archives"
        echo "  $0 extract <name> - Extract specific archive"
        echo "  $0 extract-all    - Extract all archives"
        echo "  $0 status         - Show current status"
        echo "  $0 help           - Show this help"
        echo ""
        echo "Examples:"
        echo "  $0 list"
        echo "  $0 extract L8_10_12_complete_results.tar.gz"
        echo "  $0 extract-all"
        ;;
esac

echo ""
