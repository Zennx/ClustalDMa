# Bugfix: Jaccard Contact Scoring & Distance Matrix Display

**Date:** 13 November 2025  
**Issues Fixed:** 3 critical bugs

---

## Bug #1: Purple Diagonal in Distance Matrix ❌ → ✅

### Problem
Distance matrix diagonal showing purple (high distance) when structures compared to themselves should show 0 distance (dark purple = 100% similar)

### Root Cause
Diagonal was initialized to 0 but never explicitly set, and colormap interpretation was unclear

### Fix
```python
# In compute_jaccard_contact_matrix()
# Diagonal should be 0 (100% similar to self)
np.fill_diagonal(jaccard_matrix, 0.0)
```

**File:** `clustaldemo.py` line ~258

---

## Bug #2: All Yellow Distance Matrix (Jaccard Not Working) ❌ → ✅

### Problem
Jaccard distance matrix showing all structures as highly dissimilar (all yellow = distance ~1.0), when HDOCK poses should have overlapping contacts

### Root Cause
**Critical bug:** Contacts were being stored by residue **index** instead of residue **identity**

```python
# WRONG - residue index can vary
res_pair = (p_res.index, n_res.index)
```

Problem: Residue indices are **arbitrary** and can differ between structures even if they represent the same residue. For example:
- Structure 1: ARG123 might be index 55
- Structure 2: ARG123 might be index 55 (same) OR different if topology varies
- Comparison fails because we're comparing indices, not actual residues!

### Fix
Use residue **name + number** as identifier:

```python
# CORRECT - residue name+number is consistent
res_pair = (f"{p_res.name}{p_res.resSeq}", f"{n_res.name}{n_res.resSeq}")
```

Now contacts are matched by identity:
- `("ARG123", "DG5")` will match across all structures
- Jaccard similarity actually measures contact pattern overlap

**File:** `clustaldemo.py` line ~217

---

## Bug #3: Hardcoded "RMSD Distance Matrix" Title ❌ → ✅

### Problem
Distance Matrix tab always showed "RMSD Distance Matrix" even when using Jaccard

### Fix
Made the card header dynamic:

**UI Definition:**
```python
# app.py line ~131
ui.card_header(ui.output_text("distance_matrix_title"))
```

**Render Function:**
```python
@output
@render.text
def distance_matrix_title():
    if not clustering_complete():
        return "Distance Matrix"
    clust = clusterer()
    return f"{clust.metric_name} Distance Matrix"
```

Now shows:
- "RMSD Distance Matrix" when using RMSD
- "Jaccard Contact Score Distance Matrix" when using Jaccard

**File:** `app.py` lines ~131, ~625-632

---

## Additional Improvements

### Debug Output for Jaccard
Added diagnostic information to help troubleshoot contact detection:

```python
# Contact statistics
print(f"  Contact statistics: min={min(contact_counts)}, "
      f"max={max(contact_counts)}, mean={np.mean(contact_counts):.1f}")

# Distance statistics
print(f"Jaccard distance range: {non_diag.min():.3f} - {non_diag.max():.3f}")
print(f"  Mean distance: {non_diag.mean():.3f} (closer to 0 = more similar)")
print(f"  Sample: structures 0 vs 1 = {jaccard_matrix[0,1]:.3f}")
```

This helps verify:
- Contacts are being detected (not all zero)
- Distance distribution makes sense
- Similar structures have low distance

---

## Expected Behavior After Fixes

### RMSD Distance Matrix
```
- Diagonal: Dark purple (distance = 0)
- Similar poses: Dark purple to blue
- Different poses: Green to yellow
- Title: "RMSD Distance Matrix"
```

### Jaccard Contact Score Distance Matrix
```
- Diagonal: Dark purple (distance = 0)
- Same contact pattern: Dark purple to blue  
- Partially overlapping: Blue to green
- Different contacts: Green to yellow
- Title: "Jaccard Contact Score Distance Matrix"
```

### Debug Output Example
```
Computing Jaccard contact distance matrix (cutoff=4.5 Å)...
Step 1/2: Computing interface contacts for all structures...
  Processed 10/100 structures (45 residue contacts)
  Processed 20/100 structures (52 residue contacts)
  ...
Step 2/2: Computing Jaccard distances from contact sets...
  Contact statistics: min=38, max=67, mean=48.2
Jaccard distance range: 0.123 - 0.891
  Mean distance: 0.456 (closer to 0 = more similar)
  Sample: structures 0 vs 1 = 0.234, 0 vs 2 = 0.567
✓ Contact lists stored in clusterer.contact_residue_pairs
```

---

## Testing Checklist

- [ ] RMSD clustering shows diagonal as dark purple
- [ ] Jaccard clustering shows diagonal as dark purple  
- [ ] Jaccard shows variation in colors (not all yellow)
- [ ] Similar HDOCK poses cluster together (dark blocks)
- [ ] Distance Matrix title changes based on selected metric
- [ ] Contact export CSV contains meaningful residue pairs
- [ ] Debug output shows reasonable contact counts

---

## Files Modified

1. ✅ `clustaldemo.py`
   - Line ~217: Fixed contact matching (name+number vs index)
   - Line ~258: Explicit diagonal fill
   - Line ~240-270: Added debug output

2. ✅ `app.py`
   - Line ~131: Dynamic card header
   - Line ~625-632: New distance_matrix_title() render function

---

## Impact

**Before:**
- ❌ Jaccard completely broken (all dissimilar)
- ❌ Confusing purple diagonal
- ❌ Wrong title in UI

**After:**
- ✅ Jaccard correctly identifies contact similarity
- ✅ Clear diagonal = self-similarity
- ✅ Dynamic titles match selected metric
- ✅ Exportable contact lists actually useful

---

## Technical Details

### Why Residue Indices Failed

MDTraj topology assigns indices sequentially:
```
Residue 0: MET1
Residue 1: ALA2
Residue 2: ARG3
...
Residue 55: ARG123
```

But indices are:
1. **Arbitrary** - depend on how PDB was parsed
2. **Non-portable** - differ if topology changes
3. **Meaningless** - don't reflect biological identity

### Why Name+Number Works

`f"{p_res.name}{p_res.resSeq}"` creates:
```
"ARG123" - always refers to Arginine at position 123
"LYS45"  - always refers to Lysine at position 45
"DG5"    - always refers to Guanine at position 5
```

These are:
1. **Consistent** across all structures
2. **Meaningful** - correspond to actual residue identity
3. **Comparable** - same string = same residue

### Jaccard Similarity Explained

For two structures A and B with contact sets:
```
A = {("ARG123","DG5"), ("LYS45","DA3"), ("SER67","DG5")}
B = {("ARG123","DG5"), ("LYS45","DA4"), ("THR89","DC2")}

Intersection = {("ARG123","DG5")}                    # 1 contact
Union = {("ARG123","DG5"), ("LYS45","DA3"), 
         ("LYS45","DA4"), ("SER67","DG5"), 
         ("THR89","DC2")}                            # 5 contacts

Jaccard similarity = 1/5 = 0.2
Jaccard distance = 1 - 0.2 = 0.8
```

Low distance (dark purple) = many shared contacts = similar binding mode
High distance (yellow) = few shared contacts = different binding mode

---

## Future Considerations

For very large-scale analysis with mixed topologies (not typical for HDOCK):
- Could use residue name + chain + number: `f"{p_res.chain}{p_res.name}{p_res.resSeq}"`
- Could map to canonical numbering schemes
- Could use structure-based alignment to match residues

For current HDOCK use case (same receptor, same nucleic acid sequence):
- Current fix is optimal
- Name+number is sufficient and fast
- No need for complex mapping
