"""Compatibility entry point for the automated CQUAD4 DMIG comparison.

The original version of this script duplicated geometries with mock objects
and wrote draft CSV files beside the Nastran references.  Validation now reads
the real decks for cases 51-54 and never modifies the reference directory.

New code should invoke ``run_cquad4_nastran_comparison.py`` directly.  This
wrapper is kept so existing local commands continue to work.
"""

from run_cquad4_nastran_comparison import main


if __name__ == "__main__":
    main()
