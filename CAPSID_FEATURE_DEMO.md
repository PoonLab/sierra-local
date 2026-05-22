# Capsid/Lenacapavir Resistance Support in Sierra-Local

## Overview

The `capsid` branch adds support for detecting CA (capsid) gene mutations and scoring lenacapavir (LEN) resistance.

## Features Added

### 1. CA Gene Processing
- **Gene coordinates**: 1186-1878 nt (Gag region), 1-231 AA
- **Minimum overlap**: 60 amino acids required
- **Reference frame**: Uses Gag start (790) instead of Pol start (2085)

### 2. isCapsidResistance Flag
New boolean field added to each mutation in `alignedGeneSequences`:
```json
{
  "position": 56,
  "text": "L56I",
  "isCapsidResistance": true,
  "primaryType": "Major"
}
```

### 3. Lenacapavir Resistance Scoring
- **Drug**: LEN (lenacapavir)
- **Drug Class**: CAI (Capsid Inhibitor)
- **Resistance mutations** (from HIVDB 9.8):
  - Major: L56I (60), I66I (60), Q67H (30), N/S/H70 (30/30/60), D/S74 (60/30)
  - Accessory: S57 (30), Y/N/K67 (20), R70 (10), T/E/S105 (10), N/C107 (10)

## Verified Working

✅ CA gene appears in `alignedGeneSequences`  
✅ CA gene appears in `drugResistance` with LEN scoring  
✅ `isCapsidResistance` flag present on all CA mutations  
✅ Method correctly identifies resistance mutations  

## Example Output

```json
{
  "alignedGeneSequences": [
    {
      "gene": {"name": "CA"},
      "firstAA": 1,
      "lastAA": 231,
      "mutations": [
        {
          "position": 56,
          "text": "L56I",
          "isCapsidResistance": true,
          "primaryType": "Major"
        }
      ]
    }
  ],
  "drugResistance": [
    {
      "gene": {"name": "CA"},
      "drugScores": [
        {
          "drug": {"name": "LEN", "displayAbbr": "LEN"},
          "drugClass": {"name": "CAI"},
          "score": 60,
          "level": 5,
          "text": "High-Level Resistance"
        }
      ]
    }
  ]
}
```

## Testing

Unit test verification:
```python
>>> writer.is_capsid_resistance('CA', 56, 'I')  # L56I Major
True
>>> writer.is_capsid_resistance('CA', 67, 'H')  # Q67H Major
True
>>> writer.is_capsid_resistance('CA', 105, 'E') # T105E Accessory
True
>>> writer.is_capsid_resistance('CA', 145, 'E') # Non-resistance
False
```

## Files Modified

1. `sierralocal/jsonwriter.py`
   - Added `is_capsid_resistance()` method
   - Added `isCapsidResistance` field to mutation output

2. `sierralocal/nucaminohook.py`
   - Added CA gene coordinates (1186-1878)
   - Added gag_start reference point (790)
   - Updated gene_map to handle Gag vs Pol reference frames
   - Added CA to min_overlap dictionary

## Next Steps

For Cascade pipeline integration:
1. Extract capsid resistance mutations from sierra-local output
2. Track across abundance thresholds (10%, 20%, Majority)
3. Create MultiQC table for capsid mutations
4. Add to HIV antiviral resistance reporting
