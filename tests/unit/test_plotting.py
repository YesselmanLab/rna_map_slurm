"""Tests for plotting functions."""

from __future__ import annotations

from rna_map_slurm.plotting.pop_avg import colors_for_sequence


class TestColorsForSequence:
    """Tests for colors_for_sequence function."""

    def test_adenine_red(self) -> None:
        """Test that A is red."""
        colors = colors_for_sequence("A")
        assert colors == ["red"]

    def test_cytosine_blue(self) -> None:
        """Test that C is blue."""
        colors = colors_for_sequence("C")
        assert colors == ["blue"]

    def test_guanine_orange(self) -> None:
        """Test that G is orange."""
        colors = colors_for_sequence("G")
        assert colors == ["orange"]

    def test_thymine_green(self) -> None:
        """Test that T is green."""
        colors = colors_for_sequence("T")
        assert colors == ["green"]

    def test_uracil_green(self) -> None:
        """Test that U is green."""
        colors = colors_for_sequence("U")
        assert colors == ["green"]

    def test_full_sequence(self) -> None:
        """Test a full sequence."""
        colors = colors_for_sequence("ACGU")
        assert colors == ["red", "blue", "orange", "green"]

    def test_unknown_base(self) -> None:
        """Test unknown base returns gray."""
        colors = colors_for_sequence("N")
        assert colors == ["gray"]
