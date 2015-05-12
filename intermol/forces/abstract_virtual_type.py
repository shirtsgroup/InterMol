from intermol.forces.abstract_type import AbstractType


class AbstractVirtualType(AbstractType):
    __slots__ = ['placeholder']

    virtualsitelist = ['two',
                       'three_linear',
                       'three_fd',
                       'three_fad',
                       'three_out',
                       'four_fdn',
                       'n_cog',
                       'n_com',
                       'n_cow'
                       ]

    def __init__(self, placeholder=False):
        super(AbstractVirtualType, self).__init__()
        self.placeholder = placeholder
