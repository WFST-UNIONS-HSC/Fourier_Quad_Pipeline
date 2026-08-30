#ifndef CATALOG_COMBINER_HPP
#define CATALOG_COMBINER_HPP

#include <vector>
#include <string>

namespace CatalogCombiner {
    // ==========================================
    // Function: Combine one exposure into its final source catalog
    // Method: Preserve the 1-based exposure identity ahead of the CCD column.
    // ==========================================
    void combineExpoCatalog(int nchip,
                            const std::vector<std::string>& imageFiles,
                            const std::string& dirOutput,
                            int expo_index,
                            float chi2);
    void procComb(int iexpo);
}

#endif // CATALOG_COMBINER_HPP
