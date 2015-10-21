#pragma once

namespace lepID {
    enum lepID { L, T, Count };
}

namespace lepIso {
    enum lepIso {no, L, T, Count };
}

namespace btagWP {
    enum btagWP { no, L, M, T, Count };
}

namespace jetPair {
//    enum jetPair { mh, pt, ht, btag, ptOverM, mc, Count }; //FIXME: ultimately
    enum jetPair { ht, mh, pt, csv, jp, ptOverM, Count };
}
