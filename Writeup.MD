

**Revise the mass**

    float mass_revised = mass * (1.2114 + 1.2092 + 1.2140 + 1.2140) /4; // Increase the mass by 21.215%  
    cmd.desiredThrustsN[0] = mass_revised * 9.81f / 4.f;                // * 1.2114; // front left  
    cmd.desiredThrustsN[1] = mass_revised * 9.81f / 4.f;                // * 1.2092; // front right  
    cmd.desiredThrustsN[2] = mass_revised * 9.81f / 4.f;                // * 1.2140; // rear left  
    cmd.desiredThrustsN[3] = mass_revised * 9.81f / 4.f;                // * 1.2140; // rear right  
