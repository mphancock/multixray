import IMP
import IMP.atom


class WeightThermostat(IMP.OptimizerState):
    def __init__(
            self,
            msmc_m,
            rset_xray,
            md,
            T_target,
            warmup_steps=50
    ):
        IMP.OptimizerState.__init__(self, msmc_m.get_m(), "WeightThermostat%1%")
        self.msmc_m = msmc_m
        self.m = msmc_m.get_m()
        self.rset_xray = rset_xray
        self.md = md
        self.T_target = T_target
        self.warmup_steps = warmup_steps

        self.on = True

    def get_msmc_m(self):
        return self.msmc_m

    def get_rset_xray(self):
        return self.rset_xray

    def get_md(self):
        return self.md

    def get_T_target(self):
        return self.T_target

    def get_on(self):
        return self.on

    def get_warmup_steps(self):
        return self.warmup_steps

    def set_msmc_m(
        self,
        msmc_m
    ):
        self.msmc_m = msmc_m

    def set_rset_xray(
        self,
        rset_xray
    ):
        self.rset_xray = rset_xray

    def set_md(
        self,
        md
    ):
        self.md = md

    def set_T_target(
        self,
        T_target
    ):
        self.T_target = T_target

    def set_warmup_steps(
        self,
        warmup_steps
    ):
        self.warmup_steps = warmup_steps

    def set_on(
        self,
        on
    ):
        self.on = on

    def do_update(
        self,
        call
    ):
        ## get number of updates is the number of times that do update has been called
        if self.get_on() and self.get_number_of_updates()*self.get_period() >= self.warmup_steps:
            w_xray_cur = self.rset_xray.get_restraint(0).get_weight()

            energy = self.md.get_kinetic_energy()
            T_cur = self.md.get_kinetic_temperature(ekinetic=energy)

            w_xray = w_xray_cur * self.T_target / T_cur

            for i in range(self.rset_xray.get_number_of_restraints()):
                self.rset_xray.get_restraint(i).set_weight(w_xray)

