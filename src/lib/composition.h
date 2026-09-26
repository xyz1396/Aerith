#pragma once
#include <array>
#include <cstddef>

namespace sipros
{

constexpr std::size_t ElementCount = 6;
constexpr std::size_t IsotopeSourceCount = 3;

enum class Element : std::size_t
{
	Carbon = 0,
	Hydrogen = 1,
	Oxygen = 2,
	Nitrogen = 3,
	Phosphorus = 4,
	Sulfur = 5
};

// Atom provenance is part of the compiled chemistry. Only
// Biosynthetic atoms follow the selected SIP abundance. ReagentNatural and
// DigestionSolvent atoms retain the natural isotope distribution of the real
// element when isotope envelopes are convolved.
enum class IsotopeSource : std::size_t
{
	Biosynthetic = 0,
	ReagentNatural = 1,
	DigestionSolvent = 2
};

using AtomCounts = std::array<int, ElementCount>;

struct SourcedComposition
{
	std::array<AtomCounts, IsotopeSourceCount> atoms{};

	AtomCounts &operator[](IsotopeSource source)
	{
		return atoms[static_cast<std::size_t>(source)];
	}

	const AtomCounts &operator[](IsotopeSource source) const
	{
		return atoms[static_cast<std::size_t>(source)];
	}

	SourcedComposition &operator+=(const SourcedComposition &other)
	{
		for (std::size_t source = 0; source < IsotopeSourceCount; ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				atoms[source][element] += other.atoms[source][element];
			}
		}
		return *this;
	}

	SourcedComposition &operator-=(const SourcedComposition &other)
	{
		for (std::size_t source = 0; source < IsotopeSourceCount; ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				atoms[source][element] -= other.atoms[source][element];
			}
		}
		return *this;
	}

	AtomCounts total() const
	{
		AtomCounts result{};
		for (std::size_t source = 0; source < IsotopeSourceCount; ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				result[element] += atoms[source][element];
			}
		}
		return result;
	}

	AtomCounts naturalSourceTotal() const
	{
		AtomCounts result{};
		for (std::size_t source =
				 static_cast<std::size_t>(IsotopeSource::ReagentNatural);
			 source < IsotopeSourceCount;
			 ++source)
		{
			for (std::size_t element = 0; element < ElementCount; ++element)
			{
				result[element] += atoms[source][element];
			}
		}
		return result;
	}
};

inline SourcedComposition operator+(SourcedComposition left,
									const SourcedComposition &right)
{
	left += right;
	return left;
}

inline SourcedComposition operator-(SourcedComposition left,
									const SourcedComposition &right)
{
	left -= right;
	return left;
}

inline SourcedComposition compositionFrom(IsotopeSource source,
										const AtomCounts &counts)
{
	SourcedComposition composition;
	composition[source] = counts;
	return composition;
}

} // namespace sipros


using AtomCounts = sipros::AtomCounts;
using Composition = sipros::SourcedComposition;
using IsotopeSource = sipros::IsotopeSource;
