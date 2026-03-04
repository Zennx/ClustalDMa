# Upload Workflow Fix & Development Mode Setup

## Date: 25 November 2025

## Changes Summary

### 1. Development Docker Setup ✅

Created `docker-compose.dev.yml` with:
- Mounts entire working directory to `/app`
- Auto-reload enabled via `--reload` flag
- Instant code changes without rebuilding
- Same ports: 8000 (app), 8080 (file browser)

**Launch script**: `./launch_dev.sh`

### 2. Fixed Upload/Parse Workflow ✅

**Problem**: 
- App was auto-parsing on upload
- Showed "0 PDB structures generated" message immediately
- Analysis failed because parsing happened during upload

**Solution**:
Added "Generate Poses" button to separate upload from parsing:

1. **Upload Stage**: User uploads .out files
   - Files stored in temp directory
   - Status: "✓ N .out file(s) ready - Click 'Generate Poses'"

2. **Parse Stage**: User clicks "🔄 Generate Poses"
   - HDOCK parser runs on uploaded files
   - Creates PDB structures (controlled by slider)
   - Status: "✓ N PDB structure(s) ready for analysis"

3. **Analysis Stage**: User clicks "▶️ Start Analysis"
   - Runs clustering on generated PDBs
   - No more "no PDBs found" error

### 3. Code Changes

**app_main.py**:
- Added `ui.input_action_button("generate_poses", "🔄 Generate Poses")`
- Split `handle_upload()` into two functions:
  - `handle_upload()`: Only stores .out files
  - `handle_generate_poses()`: Parses files and generates PDBs
- Updated `upload_status()`:
  - Shows "ready - Click 'Generate Poses'" when files uploaded
  - Shows "ready for analysis" when PDBs generated

**New Files**:
- `docker-compose.dev.yml`: Development mode configuration
- `launch_dev.sh`: Quick launch script for dev mode
- `DEV_WORKFLOW.md`: Documentation for development workflow

## Usage

### Development Mode (Current/Testing)
```bash
./launch_dev.sh
# Code changes auto-reload - just refresh browser!
```

### Production Mode (When Ready)
```bash
./launch_docker.sh
# Rebuild required for code changes
```

## Workflow Now

1. Upload .out files → Status shows "ready - Click 'Generate Poses'"
2. Click "🔄 Generate Poses" → Parsing happens, PDbs created
3. Click "▶️ Start Analysis" → Clustering runs
4. Click "📥 Download Results" → Get ZIP with all outputs

## Benefits

✅ No more tedious rebuilds during development
✅ Clear separation: Upload → Parse → Analyze
✅ No confusing "0 structures" messages
✅ User controls when parsing happens
✅ Fast iteration for testing changes

## Testing Status

- [x] Development mode launches successfully
- [x] Generate Poses button appears in UI
- [x] Upload handler only stores files (no parsing)
- [x] Status messages updated correctly
- [ ] Test actual .out file upload
- [ ] Test Generate Poses button functionality
- [ ] Test full workflow end-to-end
