# Reaction Contribution Web App
This small static web app allows you to load up a reactions contributions JSON file and visualize the contributions of different reactions to different species. This is useful for debugging.

The JSON file must have the schema

```json
{
  i: {
     t: float,
     dt: float,
     reaction_contribution: {
        Species: {
           reaction: float
        }
     }
  }
}

```
where i counts from 0 -> the number of timesteps, Species is the species name (so if at one timestep there are 10 species tracked there will be 10 species entries), and reaction is the projectile-ejectile name of a reaction. 
