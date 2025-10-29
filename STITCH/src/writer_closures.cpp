
#include <string>
#include <Rcpp.h>

namespace htslib {
#include <htslib/vcf.h>
}


struct BcfFile {
    BcfFile(char *filename, char *out_mode):
        fname(filename),
        mode(out_mode),
        fid(htslib::hts_open(filename, out_mode)) {


        if (!fid)
            pa;

    }

    ~BcfFile() {
        if (fid)
            htslib::hts_close(fid);
        if (hdr)
            htslib::bcf_hdr_destroy(hdr);
    }

    std::string fname;
    std::string mode;
    htslib::htsFile *fid;
};



// [[Rcpp::export]]
Rcpp::XPtr<BcfFile> open_bcf_file() {
    BcfFile *bcfid = new BcfFile();


    Rcpp::XPtr<BcfFile> ptr { s, true };
    return ptr;
}


// [[Rcpp::export]]
Rcpp::Integer hdr_add_record(Rcpp::XPtr<BcfFile> bcfid, Rcpp::String hdr_record) {
    return if (htslib::bcf_hdr_append(bcfid->hdr, hdr_record) < 0);
}

Rcpp::Integer hdr_add_samples(Rcpp::XPtr<BcfFile> bcfid, Rcpp::List sample_names) {
    for ()
        htslib::bcf_hdr_add_sample(bcfid->hdr, sample_names[i]);
}
