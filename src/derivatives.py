import time

import IMP
import IMP.core
import IMP.algebra


## function to get the derivatives of a restraint r with respect to the position of pids
def evaluate_df_dict(
        m,
        pids,
        r
):
    ## a call to evaluate will override the current derivative values and the derivatives must be reset to their original values afterwards
    ## dictionary to store restraint derivatives
    df_dict = dict()

    ## dictionary to store original derivatives
    df_0_dict = dict()

    ## store the original derivatives
    for pid in pids:
        xyz = IMP.core.XYZ(m, pid)
        df_0 = xyz.get_derivatives()
        df_0_dict[pid] = df_0

    # don't need to 0 because calling the top-level evaluate method which sets derivatives to 0 first
    r.evaluate(True)

    ## get the restraint derivatives and reset to the original derivatives
    for pid in pids:
        xyz = IMP.core.XYZ(m, pid)
        df_dict[pid] = xyz.get_derivatives()
        set_df(m, pid, df_0_dict[pid])

    return df_dict


# get_df_dict returns a dictionary with a particle index as the key and the corresponding 3D first derivative as the entry.
def get_df_dict(
        m,
        pids,
        r
):
    # print("r: {}".format(r.get_name()))
    # See if the restraint has a get_df call that stores the 3d vector first derivatives of all atoms.
    try:
        # get_df_dict returns stored a dictionary where keys are particle indexes and entries are the corresponding derivatives (a 3D vector).
        df_dict = r.get_df_dict()
    except AttributeError:
        # Otherwise, the derivative dictionary must be constructed.
        df_dict = evaluate_df_dict(
            m=m,
            pids=pids,
            r=r
        )

    return df_dict

def set_df(
        m,
        pid,
        df
):
    da = IMP.DerivativeAccumulator(1)
    xyz = IMP.core.XYZ(m, pid)
    df_0 = xyz.get_derivatives()
    xyz.add_to_derivatives(
        v=df-df_0,
        d=da
    )


# r_1 and r_2 must be IMP restraints or restraint sets.
def get_df_mag_ratio(
        m,
        pids,
        r1,
        r2
):
    df_dict_1 = r1.get_df_dict()
    df_dict_2 = r2.get_df_dict()

    # dx_dict_1 = get_df_dict(
    #     m=m,
    #     pids=pids,
    #     r=r1
    # )
    # dx_dict_2 = get_df_dict(
    #     m=m,
    #     pids=pids,
    #     r=r2
    # )
    avg_mag_1, avg_mag_2 = 0, 0

    for pid in pids:
        df_1 = df_dict_1[pid]
        df_2 = df_dict_2[pid]

        avg_mag_1 = avg_mag_1 + df_1.get_magnitude()
        avg_mag_2 = avg_mag_2 + df_2.get_magnitude()

    ## divide by the sums
    mag_ratio = avg_mag_1 / avg_mag_2

    return mag_ratio


# The optimizer state will automatically adjust the weight of restraint 2 to match restraint 1.
class RestraintWeightOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            pids,
            r1,
            r2
    ):
        IMP.OptimizerState.__init__(self, m, "RestraintWeightOptimizerState%1%")
        self.m = m
        self.pids = pids
        self.r1 = r1
        self.r2 = r2

    def do_update(self, call):
        mag_ratio = get_df_mag_ratio(
            self.m,
            self.pids,
            self.r1,
            self.r2
        )

        # self.r2.set_weight(mag_ratio)
        print("weight: {}".format(self.r2.get_weight()))


