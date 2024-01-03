<template>
  <PageHeader highlighted-word="Molecule" other-words="Classifier"/>
  <ChemDrawTool @smiles="storeSmiles"/>
  <v-btn :color="$route.meta.navBarColor" v-if="smiles" @click="displayClassifications">Get classifications</v-btn>
  <v-chip-group v-if="classifications">
    <v-chip v-for="classification in classifications" :key="classification">
      {{ classification }}
    </v-chip>
  </v-chip-group>
</template>

<script>
import axios from 'axios'
import ChemDrawTool from '@/components/ChemDrawTool.vue'
import PageHeader from '@/components/PageHeader.vue'

export default {
  components: { ChemDrawTool, PageHeader },
  data () {
    return {
      smiles: null,
      classifications: null
    }
  },
  methods: {
    storeSmiles (smiles) {
      this.smiles = smiles
      this.classifications = null
    },
    async displayClassifications () {
      await axios.get('/mol/classifications', { params: { mol_smiles: this.smiles } })
        .then(res => {
          this.classifications = res.data
        })
        .catch(error => {
          console.log(error.response)
        })
    }
  }
}
</script>
