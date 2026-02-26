-- ============================================================================
-- Migration: 004_add_chainname_classifications.sql
-- Description: Add chain name label columns based on chain ID
-- 
-- This migration adds protein name labels for each chain ID:
-- 1. chainID_label_i: protein name for chainID_i
-- 2. chainID_label_j: protein name for chainID_j
-- 
-- Mappings are derived from chain_labels.yaml
-- ============================================================================

-- ----------------------------------------------------------------------------
-- SECTION 1: Chain Name Label (chainID_label_i)
-- Purpose: Map first chain ID to its protein name
-- ----------------------------------------------------------------------------

ALTER TABLE data ADD COLUMN chainID_label_i TEXT;

UPDATE data
SET chainID_label_i = CASE chainID_i
    -- CP24
    WHEN '8' THEN 'CP24'
    WHEN '4' THEN 'CP24'

    -- CP26
    WHEN 's' THEN 'CP26'

    -- CP29
    WHEN 'r' THEN 'CP29'
    WHEN 'R' THEN 'CP29'

    -- CP43
    WHEN 'c' THEN 'CP43'

    -- LHCBM
    WHEN '5' THEN 'LHCBM'
    WHEN '6' THEN 'LHCBM'
    WHEN 'g' THEN 'LHCBM'
    WHEN 'n' THEN 'LHCBM'
    WHEN 'Y' THEN 'LHCBM'

    -- LHCB3
    WHEN '7' THEN 'LHCB3'

    -- PsbS
    WHEN '9' THEN 'PsbS'
    WHEN 'A' THEN 'PsbS'

    -- PsbA
    WHEN 'a' THEN 'PsbA'

    -- PsbB
    WHEN 'B' THEN 'PsbB'
    WHEN 'b' THEN 'PsbB'

    -- PsbE
    WHEN 'E' THEN 'PsbE'
    WHEN 'e' THEN 'PsbE'

    -- PsbF
    WHEN 'F' THEN 'PsbF'
    WHEN 'f' THEN 'PsbF'

    -- PsbG
    WHEN 'G' THEN 'PsbG'

    -- PsbH
    WHEN 'H' THEN 'PsbH'
    WHEN 'h' THEN 'PsbH'

    -- PsbI
    WHEN 'I' THEN 'PsbI'
    WHEN 'i' THEN 'PsbI'

    -- PsbJ
    WHEN 'J' THEN 'PsbJ'
    WHEN 'j' THEN 'PsbJ'

    -- PsbK
    WHEN 'K' THEN 'PsbK'
    WHEN 'k' THEN 'PsbK'

    -- PsbL
    WHEN 'L' THEN 'PsbL'
    WHEN 'l' THEN 'PsbL'

    -- PsbM
    WHEN 'M' THEN 'PsbM'
    WHEN 'm' THEN 'PsbM'

    -- PsbO
    WHEN 'o' THEN 'PsbO'

    -- PsbP
    WHEN 'P' THEN 'PsbP'
    WHEN 'p' THEN 'PsbP'

    -- PsbQ
    WHEN 'Q' THEN 'PsbQ'
    WHEN 'q' THEN 'PsbQ'

    -- PsbT
    WHEN 'T' THEN 'PsbT'
    WHEN 't' THEN 'PsbT'

    -- PsbW
    WHEN 'W' THEN 'PsbW'
    WHEN 'w' THEN 'PsbW'

    -- PsbX
    WHEN 'X' THEN 'PsbX'
    WHEN 'x' THEN 'PsbX'

    -- PsbZ
    WHEN 'Z' THEN 'PsbZ'
    WHEN 'z' THEN 'PsbZ'

    ELSE 'unknown'
END;


-- ----------------------------------------------------------------------------
-- SECTION 2: Chain Name Label (chainID_label_j)
-- Purpose: Map second chain ID to its protein name
-- ----------------------------------------------------------------------------

ALTER TABLE data ADD COLUMN chainID_label_j TEXT;

UPDATE data
SET chainID_label_j = CASE chainID_j
    -- CP24
    WHEN '8' THEN 'CP24'
    WHEN '4' THEN 'CP24'

    -- CP26
    WHEN 's' THEN 'CP26'

    -- CP29
    WHEN 'r' THEN 'CP29'
    WHEN 'R' THEN 'CP29'

    -- CP43
    WHEN 'c' THEN 'CP43'

    -- LHCBM
    WHEN '5' THEN 'LHCBM'
    WHEN '6' THEN 'LHCBM'
    WHEN 'g' THEN 'LHCBM'
    WHEN 'n' THEN 'LHCBM'
    WHEN 'Y' THEN 'LHCBM'

    -- LHCB3
    WHEN '7' THEN 'LHCB3'

    -- PsbS
    WHEN '9' THEN 'PsbS'
    WHEN 'A' THEN 'PsbS'

    -- PsbA
    WHEN 'a' THEN 'PsbA'

    -- PsbB
    WHEN 'B' THEN 'PsbB'
    WHEN 'b' THEN 'PsbB'

    -- PsbE
    WHEN 'E' THEN 'PsbE'
    WHEN 'e' THEN 'PsbE'

    -- PsbF
    WHEN 'F' THEN 'PsbF'
    WHEN 'f' THEN 'PsbF'

    -- PsbG
    WHEN 'G' THEN 'PsbG'

    -- PsbH
    WHEN 'H' THEN 'PsbH'
    WHEN 'h' THEN 'PsbH'

    -- PsbI
    WHEN 'I' THEN 'PsbI'
    WHEN 'i' THEN 'PsbI'

    -- PsbJ
    WHEN 'J' THEN 'PsbJ'
    WHEN 'j' THEN 'PsbJ'

    -- PsbK
    WHEN 'K' THEN 'PsbK'
    WHEN 'k' THEN 'PsbK'

    -- PsbL
    WHEN 'L' THEN 'PsbL'
    WHEN 'l' THEN 'PsbL'

    -- PsbM
    WHEN 'M' THEN 'PsbM'
    WHEN 'm' THEN 'PsbM'

    -- PsbO
    WHEN 'o' THEN 'PsbO'

    -- PsbP
    WHEN 'P' THEN 'PsbP'
    WHEN 'p' THEN 'PsbP'

    -- PsbQ
    WHEN 'Q' THEN 'PsbQ'
    WHEN 'q' THEN 'PsbQ'

    -- PsbT
    WHEN 'T' THEN 'PsbT'
    WHEN 't' THEN 'PsbT'

    -- PsbW
    WHEN 'W' THEN 'PsbW'
    WHEN 'w' THEN 'PsbW'

    -- PsbX
    WHEN 'X' THEN 'PsbX'
    WHEN 'x' THEN 'PsbX'

    -- PsbZ
    WHEN 'Z' THEN 'PsbZ'
    WHEN 'z' THEN 'PsbZ'

    ELSE 'unknown'
END;

-- ============================================================================
-- End of migration 004
-- ============================================================================
