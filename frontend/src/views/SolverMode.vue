<template>
  <PageHeader highlighted-word="Solver" other-words="Mode"/>
  <ChemDrawTool @smiles="storeSmiles"/>

</template>

<script>

import ChemDrawTool from '@/components/ChemDrawTool.vue'
import PageHeader from '@/components/PageHeader.vue'

export default {
  name: 'SolverMode',
  components: { ChemDrawTool, PageHeader },
  data () {
    return {
      smiles: null
    }
  },
  methods: {
    storeSmiles (smiles) {
      this.smiles = smiles
    }
  }

}
</script>
