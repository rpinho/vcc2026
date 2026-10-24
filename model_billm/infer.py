"""Predict source-only directions from the trained three-seed ensemble."""
import common as rt
import numpy as np,torch
from dataset import Sources
from network import ResponseMLP,blend

class Ensemble:
    def __init__(self,fold='production',device='cpu'):
        self.src=Sources([] if fold=='production' else [fold]);self.device=device;self.models=[]
        for seed in [1019,2027,4093]:
            p=rt.ROOT/f'checkpoints/{fold}/seed{seed}_weights.pt'
            model=ResponseMLP(self.src.g).to(device)
            model.load_state_dict(torch.load(p,map_location=device,weights_only=True));model.eval();self.models.append(model)
    def predict(self,targets,control,coefficient=.3):
        features=np.stack([self.src.features(t,control) for t in targets]);x=torch.from_numpy(features).to(self.device)
        with torch.no_grad(): learned=torch.stack([model(x) for model in self.models]).mean(0).cpu().numpy()
        response=blend(features[:,:self.src.g],learned,coefficient)
        return dict(response=response,pooled=features[:,:self.src.g],learned=learned,genes=self.src.genes)
