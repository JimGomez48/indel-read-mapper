__author__ = 'james'


class ReadMapper:
    __subseq_map__ = {}
    __read_map__ = {}

    def __init__(self, ref):
        self.__build_subsequence_map__(ref)

    def map_reads(self, reads):
        # TODO map reads using subsequence map
        return self.__read_map__

    def __get_best_match__(self, ):
        pass

    def __build_subsequence_map__(self, ref):
        # TODO build the subsequence map for fast read mapping
        pass